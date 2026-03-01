"""
Vercel Serverless API for OmniOmics AI
"""

import json
import os
import requests
from http.server import BaseHTTPRequestHandler
from urllib.parse import parse_qs, urlparse

# Environment variables (set in Vercel project settings)
GEMINI_API_KEY = os.environ.get("GEMINI_API_KEY", "")
KILO_API_KEY = os.environ.get("KILO_API_KEY", "")


def call_gemini(prompt: str, max_tokens: int = 2048) -> str:
    """Call Gemini API"""
    url = "https://generativelanguage.googleapis.com/v1beta/models/gemini-2.0-flash:generateContent"
    payload = {
        "contents": [{"parts": [{"text": prompt}]}],
        "generationConfig": {
            "maxOutputTokens": max_tokens,
            "temperature": 0.7,
        },
    }
    response = requests.post(
        url, params={"key": GEMINI_API_KEY}, json=payload, timeout=30
    )
    if response.status_code == 200:
        result = response.json()
        return result["candidates"][0]["content"]["parts"][0]["text"]
    return f"Gemini Error: {response.status_code} - {response.text[:200]}"


def call_kilo(prompt: str, max_tokens: int = 2048) -> str:
    """Call Kilo API (Claude)"""
    url = "https://api.kilo.ai/api/gateway/chat/completions"
    headers = {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {KILO_API_KEY}",
    }
    payload = {
        "model": "anthropic/claude-sonnet-4.5",
        "messages": [{"role": "user", "content": prompt}],
        "max_tokens": max_tokens,
        "temperature": 0.7,
    }
    response = requests.post(url, headers=headers, json=payload, timeout=30)
    if response.status_code == 200:
        result = response.json()
        return result["choices"][0]["message"]["content"]
    return f"Kilo Error: {response.status_code} - {response.text[:200]}"


def generate_analysis(prompt: str, provider: str = "gemini") -> str:
    """Generate analysis using specified provider"""
    if provider == "kilo":
        if not KILO_API_KEY:
            return (
                "KILO_API_KEY not configured. Please set it in Vercel project settings."
            )
        return call_kilo(prompt)
    else:
        if not GEMINI_API_KEY:
            return "GEMINI_API_KEY not configured. Please set it in Vercel project settings."
        return call_gemini(prompt)


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
            response = {
                "status": "ok",
                "service": "OmniOmics AI",
                "gemini": "configured" if GEMINI_API_KEY else "missing",
                "kilo": "configured" if KILO_API_KEY else "missing",
            }
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
            prompt = data.get("prompt", "")
            provider = data.get("provider", "gemini")
            file_name = data.get("fileName", "unknown")

            if not prompt:
                response = {"error": "No prompt provided"}
            else:
                # Check if it's a proteomics file analysis
                if "proteomics" in prompt.lower() or "protein" in prompt.lower():
                    analysis_prompt = f"""You are an expert bioinformatics analyst specializing in proteomics.

Analyze the following proteomics data and provide:
1. Key findings and patterns
2. Potential biological significance
3. Pathway implications
4. Recommendations for follow-up

Data to analyze:
{prompt}

Provide a detailed scientific analysis."""
                else:
                    analysis_prompt = f"""You are an expert bioinformatics analyst specializing in proteomics, transcriptomics, and multi-omics integration.

{prompt}

Provide a detailed, scientifically accurate response."""

                analysis = generate_analysis(analysis_prompt, provider)

                response = {
                    "provider": provider,
                    "fileName": file_name,
                    "analysis": analysis,
                    "status": "success",
                }
        else:
            response = {"error": "Endpoint not found"}

        self.wfile.write(json.dumps(response).encode())


# For local testing
if __name__ == "__main__":
    from http.server import HTTPServer

    # Set keys for local testing
    os.environ["GEMINI_API_KEY"] = "AIzaSyDDHHT96BJ1dnEriQmR_bkZgzZuMKEI3aI"
    os.environ["KILO_API_KEY"] = (
        "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJlbnYiOiJwcm9kdWN0aW9uIiwia2lsb1VzZXJJZCI6IjdmN2MzZDA1LTM5NWQtNDhhYS05YjkzLTgzMjIxYzRlODIwMyIsImFwaVRva2VuUGVwcGVyIjoiZjEwODA0OTItZTliMi00N2QxLThlMzktYzkwOGE1MmU1NWJjIiwidmVyc2lvbiI6MywiaWF0IjoxNzcyMzIwMjgzLCJleHAiOjE5MzAwMDAyODN9.RoQt4MKVX4oiENfou8pZA9vhihtNDtQBLVAciuS_GXE"
    )

    server = HTTPServer(("localhost", 8000), handler)
    print("Server running on http://localhost:8000")
    print("API: /api/health, /api/analyze")
    server.serve_forever()
