"""
Chat agent — lightweight Q&A mode for the GRN benchmark assistant.

Uses a plain chat completion loop (OpenAI or Anthropic) with:
  - All know-how docs injected directly into the system prompt.
  - search_viash available as a function-calling tool for Viash/Docker/Nextflow
    syntax questions not covered by the know-how docs.

This is the default agent mode. Use A1 mode (--agent a1) only when you need
to execute tasks (run Viash commands, modify files, etc.).
"""

import json
import os
import re
import sys
from pathlib import Path


# ─── Viash search (same logic as agent.py) ────────────────────────────────────

_VIASH_URL_CACHE: dict = {"urls": None}

_SEARCH_VIASH_SCHEMA = {
    "type": "function",
    "function": {
        "name": "search_viash",
        "description": (
            "Search the official Viash documentation (viash.io/guide and viash.io/reference) "
            "for Viash component configuration, Docker engine setup, Nextflow runner, argument "
            "types, resources, test benches, and CLI commands. Only call this when the question "
            "is specifically about Viash/Docker/Nextflow syntax or config details not already "
            "answered by the know-how documents in your context."
        ),
        "parameters": {
            "type": "object",
            "properties": {
                "query": {"type": "string", "description": "Question about Viash, Docker, or Nextflow."}
            },
            "required": ["query"],
        },
    },
}


def _search_viash(query: str) -> str:
    import requests
    from bs4 import BeautifulSoup

    def _get_urls():
        if _VIASH_URL_CACHE["urls"] is not None:
            return _VIASH_URL_CACHE["urls"]
        try:
            r = requests.get("https://viash.io/sitemap.xml", timeout=15)
            r.raise_for_status()
            all_urls = re.findall(r"<loc>(.*?)</loc>", r.text)
            urls = [
                u for u in all_urls
                if any(x in u for x in ["/guide/", "/reference/", "/faq"])
                and "/deprecated/" not in u
            ]
            _VIASH_URL_CACHE["urls"] = urls
            return urls
        except Exception:
            return []

    def _fetch_page(session, url):
        try:
            resp = session.get(url, timeout=10)
            if not resp.ok:
                return None
            soup = BeautifulSoup(resp.text, "html.parser")
            for tag in soup.find_all(["nav", "header", "footer", "aside", "script", "style"]):
                tag.decompose()
            main = (
                soup.find("main")
                or soup.find("article")
                or soup.find("div", class_=re.compile(r"content|prose|doc"))
            )
            text = (main or soup).get_text(separator="\n", strip=True)
            return text if len(text) >= 100 else None
        except Exception:
            return None

    try:
        keywords = [w for w in re.split(r"\W+", query.lower()) if len(w) > 2]
        urls = _get_urls()
        if not urls:
            return "Could not fetch viash.io sitemap. Check your internet connection."
        scored = sorted(urls, key=lambda u: sum(1 for kw in keywords if kw in u.lower()), reverse=True)
        session = requests.Session()
        session.headers["User-Agent"] = "Mozilla/5.0 (viash-docs-searcher)"
        parts = []
        for url in scored:
            text = _fetch_page(session, url)
            if text:
                parts.append(f"[Source: {url}]\n{text[:3000]}")
            if len(parts) >= 3:
                break
        return "\n\n---\n\n".join(parts) if parts else "No relevant viash.io pages found."
    except Exception as e:
        return f"Error searching viash.io: {e}"


# ─── System prompt builder ─────────────────────────────────────────────────────

def _build_system_prompt(know_how_dir: Path) -> str:
    parts = [
        "You are an expert assistant for the `task_grn_inference` GRN benchmarking framework (GeneRNBI).\n"
        "Answer questions about GRN methods, datasets, evaluation metrics, output format, "
        "adding new methods, and running Viash commands.\n\n"
        "The know-how documents below are your primary knowledge source — always answer from them first.\n"
        "Only call `search_viash` for Viash/Docker/Nextflow syntax details not covered here.\n"
        "Never hallucinate. If something is not in the docs, say so.\n\n"
        "─────────────────────────────────────────\n"
        "KNOW-HOW DOCUMENTS\n"
        "─────────────────────────────────────────\n"
    ]
    if know_how_dir.exists():
        for md_path in sorted(know_how_dir.glob("*.md")):
            content = md_path.read_text(encoding="utf-8")
            parts.append(f"### {md_path.stem.replace('_', ' ').title()}\n\n{content}\n")
    return "\n".join(parts)


# ─── Chat agent class ──────────────────────────────────────────────────────────

class ChatAgent:
    def __init__(self, provider: str, model: str, know_how_dir: Path):
        self._provider = provider.lower()
        self._model = model
        self._system_prompt = _build_system_prompt(know_how_dir)
        self._history: list = []

    def _make_chat_client(self):
        from openai import OpenAI
        if self._provider == "openrouter":
            return OpenAI(
                base_url="https://openrouter.ai/api/v1",
                api_key=os.environ.get("OPENROUTER_API_KEY"),
            )
        return OpenAI()

    def _call_chat(self, messages: list) -> tuple[str, list | None]:
        client = self._make_chat_client()
        response = client.chat.completions.create(
            model=self._model,
            messages=messages,
            tools=[_SEARCH_VIASH_SCHEMA],
            tool_choice="auto",
        )
        msg = response.choices[0].message
        tool_calls = msg.tool_calls or []
        return msg.content or "", tool_calls

    def _call_anthropic(self, messages: list) -> tuple[str, list | None]:
        import anthropic
        client = anthropic.Anthropic()
        tool_def = {
            "name": "search_viash",
            "description": _SEARCH_VIASH_SCHEMA["function"]["description"],
            "input_schema": _SEARCH_VIASH_SCHEMA["function"]["parameters"],
        }
        response = client.messages.create(
            model=self._model,
            max_tokens=4096,
            system=self._system_prompt,
            messages=messages,
            tools=[tool_def],
        )
        text = "".join(b.text for b in response.content if hasattr(b, "text"))
        tool_uses = [b for b in response.content if b.type == "tool_use"]
        return text, tool_uses

    def ask(self, question: str) -> str:
        if self._provider in ("openai", "openrouter"):
            return self._ask_chat(question)
        return self._ask_anthropic(question)

    def _ask_chat(self, question: str) -> str:
        self._history.append({"role": "user", "content": question})
        messages = [{"role": "system", "content": self._system_prompt}] + self._history

        for _ in range(5):  # max tool-call rounds
            content, tool_calls = self._call_chat(messages)

            if not tool_calls:
                self._history.append({"role": "assistant", "content": content})
                return content

            # Execute tool calls
            messages.append({"role": "assistant", "content": content, "tool_calls": [
                {"id": tc.id, "type": "function", "function": {"name": tc.function.name, "arguments": tc.function.arguments}}
                for tc in tool_calls
            ]})
            for tc in tool_calls:
                args = json.loads(tc.function.arguments)
                print(f"  [tool] search_viash({args.get('query', '')})", flush=True)
                result = _search_viash(args["query"])
                messages.append({"role": "tool", "tool_call_id": tc.id, "content": result})

        # Final call after tool results
        content, _ = self._call_chat(messages)
        self._history.append({"role": "assistant", "content": content})
        return content

    def _ask_anthropic(self, question: str) -> str:
        self._history.append({"role": "user", "content": question})
        messages = list(self._history)

        for _ in range(5):
            text, tool_uses = self._call_anthropic(messages)

            if not tool_uses:
                self._history.append({"role": "assistant", "content": text})
                return text

            # Build assistant message with all content blocks
            import anthropic
            client = anthropic.Anthropic()
            tool_def = {
                "name": "search_viash",
                "description": _SEARCH_VIASH_SCHEMA["function"]["description"],
                "input_schema": _SEARCH_VIASH_SCHEMA["function"]["parameters"],
            }
            # Reconstruct full content for assistant turn
            full_response = client.messages.create(
                model=self._model,
                max_tokens=4096,
                system=self._system_prompt,
                messages=messages,
                tools=[tool_def],
            )
            messages.append({"role": "assistant", "content": full_response.content})

            tool_results = []
            for tu in tool_uses:
                print(f"  [tool] search_viash({tu.input.get('query', '')})", flush=True)
                result = _search_viash(tu.input["query"])
                tool_results.append({
                    "type": "tool_result",
                    "tool_use_id": tu.id,
                    "content": result,
                })
            messages.append({"role": "user", "content": tool_results})

        text, _ = self._call_anthropic(messages)
        self._history.append({"role": "assistant", "content": text})
        return text


def initialize_chat_agent(model: str = None, provider: str = None) -> "ChatAgent":
    from dotenv import load_dotenv
    from config import LLM_PROVIDER, MODEL_ID
    from rag_builder import _missing_key_error

    load_dotenv()
    prov = (provider or LLM_PROVIDER).lower()
    mod = model or MODEL_ID

    if prov == "openai" and not os.environ.get("OPENAI_API_KEY"):
        _missing_key_error("OPENAI_API_KEY")
    elif prov == "anthropic" and not os.environ.get("ANTHROPIC_API_KEY"):
        _missing_key_error("ANTHROPIC_API_KEY")
    elif prov == "openrouter" and not os.environ.get("OPENROUTER_API_KEY"):
        _missing_key_error("OPENROUTER_API_KEY")

    know_how_dir = Path(__file__).parent.parent / "know_how"
    print(f"Creating chat agent ({prov}/{mod})...")
    agent = ChatAgent(provider=prov, model=mod, know_how_dir=know_how_dir)
    print("Chat agent ready!")
    return agent
