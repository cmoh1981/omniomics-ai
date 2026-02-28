"""
Vercel Serverless API for OmniOmics AI
"""

import json
import os
from http.server import BaseHTTPRequestHandler
from urllib.parse import parse_qs, urlparse

# Initialize providers
gemini_api_key = os.environ.get("GEMINI_API_KEY", "")
kilo_api_key = os.environ.get("KILO_API_KEY", "")


class handler(BaseHTTPRequestHandler):
    """Vercel Python handler"""

    def send_cors_headers(self):
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")

    def do_OPTIONS(self):
        self.send_response(200)
        self.send_cors_headers()
        self.end_headers()

    def do_GET(self):
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_cors_headers()
        self.end_headers()

        parsed = urlparse(self.path)

        if parsed.path == "/api/health":
            response = {"status": "ok", "service": "OmniOmics AI"}
        elif parsed.path == "/api/models":
            response = {
                "providers": ["gemini", "kilo"],
                "gemini": "gemini-2.0-flash",
                "kilo": "anthropic/claude-sonnet-4.5",
            }
        else:
            response = {
                "message": "OmniOmics AI API",
                "endpoints": {
                    "POST /api/analyze": "Analyze omics data with LLM",
                    "GET /api/health": "Health check",
                    "GET /api/models": "Available LLM models",
                },
            }

        self.wfile.write(json.dumps(response).encode())

    def do_POST(self):
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_cors_headers()
        self.end_headers()

        content_length = int(self.headers.get("Content-Length", 0))
        body = self.rfile.read(content_length).decode() if content_length > 0 else "{}"

        try:
            data = json.loads(body) if body else {}
        except:
            data = {}

        parsed = urlparse(self.path)

        if parsed.path == "/api/analyze":
            # LLM Analysis endpoint
            prompt = data.get("prompt", "")
            provider = data.get("provider", "gemini")

            if not prompt:
                response = {"error": "No prompt provided"}
            else:
                # Simple analysis response
                response = {
                    "provider": provider,
                    "prompt": prompt,
                    "analysis": f"Analysis requested: {prompt[:100]}...",
                    "status": "success",
                }
        elif parsed.path == "/api/proteomics":
            response = {
                "message": "Proteomics analysis endpoint",
                "status": "coming_soon",
            }
        else:
            response = {"error": "Endpoint not found"}

        self.wfile.write(json.dumps(response).encode())


# For local testing
if __name__ == "__main__":
    from http.server import HTTPServer

    server = HTTPServer(("localhost", 8000), handler)
    print("Server running on http://localhost:8000")
    server.serve_forever()
