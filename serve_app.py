"""Local development server for ADC Input Model.

It serves the static application and accepts generated PDF reports from the
browser so the tool can write reports to the project workspace.
"""

from __future__ import annotations

import argparse
import json
import re
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any


ROOT_DIR = Path(__file__).resolve().parent
REPORTS_DIR = ROOT_DIR / "reports"
MAX_REPORT_BYTES = 50 * 1024 * 1024


def safe_report_name(value: str) -> str:
    """Return a conservative PDF file name for a browser-submitted report."""
    name = Path(value or "").name
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("._")
    if not name:
        name = "adc-input-report.pdf"
    if not name.lower().endswith(".pdf"):
        name = f"{name}.pdf"
    return name


class ADCModelHandler(SimpleHTTPRequestHandler):
    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, directory=str(ROOT_DIR), **kwargs)

    def _send_json(self, status: int, payload: dict[str, Any]) -> None:
        body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_POST(self) -> None:
        if self.path.split("?", 1)[0] != "/api/save-report":
            self._send_json(404, {"ok": False, "error": "Unknown endpoint"})
            return

        try:
            length = int(self.headers.get("Content-Length", "0"))
        except ValueError:
            self._send_json(400, {"ok": False, "error": "Invalid content length"})
            return

        if length <= 0:
            self._send_json(400, {"ok": False, "error": "Empty report"})
            return
        if length > MAX_REPORT_BYTES:
            self._send_json(413, {"ok": False, "error": "Report is too large"})
            return

        data = self.rfile.read(length)
        if not data.startswith(b"%PDF"):
            self._send_json(400, {"ok": False, "error": "Payload is not a PDF"})
            return

        REPORTS_DIR.mkdir(exist_ok=True)
        filename = safe_report_name(self.headers.get("X-Report-File", ""))
        path = REPORTS_DIR / filename
        path.write_bytes(data)

        self._send_json(
            200,
            {
                "ok": True,
                "name": filename,
                "path": path.as_posix(),
                "bytes": len(data),
            },
        )


def main() -> None:
    global REPORTS_DIR

    parser = argparse.ArgumentParser(description="Serve ADC Input Model locally")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument(
        "--reports-dir",
        default=str(REPORTS_DIR),
        help="Directory where generated PDF reports are written",
    )
    args = parser.parse_args()
    reports_dir = Path(args.reports_dir).expanduser()
    if not reports_dir.is_absolute():
        reports_dir = ROOT_DIR / reports_dir
    REPORTS_DIR = reports_dir.resolve()

    server = ThreadingHTTPServer((args.host, args.port), ADCModelHandler)
    print(f"Serving ADC Input Model at http://{args.host}:{args.port}/")
    print(f"PDF reports will be written to {REPORTS_DIR}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping server.")
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
