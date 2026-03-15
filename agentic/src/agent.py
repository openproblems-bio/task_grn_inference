"""
GeneRNBI orchestrating agent — powered by Biomni A1.

One tool (registered via a1.add_tool() — standard Biomni API):
  - search_viash:  Live page fetch from viash.io/guide + /reference (no pre-indexing)

Know-how documents in agentic/know_how/ are injected into A1's know_how_loader
so configure() includes them in the system prompt automatically.
"""

import os
import re
from pathlib import Path

from dotenv import load_dotenv

from config import LLM_PROVIDER, MODEL_ID
from rag_builder import _missing_key_error

_SOURCE_MAP = {"openai": "OpenAI", "anthropic": "Anthropic", "openrouter": "OpenAI"}

# Module-level state — set by initialize_agent() before add_tool() calls so
# that inspect.getsource() on the top-level functions works cleanly.
_VIASH_URL_CACHE: dict = {"urls": None}


# ─── Tool functions (module-level so inspect.getsource captures them correctly) ──

def search_viash(query: str) -> str:
    """Search the official Viash documentation (viash.io/guide and viash.io/reference)
    for information about Viash component configuration, Docker engine setup,
    Nextflow runner, argument types, resources, test benches, and CLI commands.
    Use this tool when the user has a Viash, Docker, or Nextflow related question
    during integration or troubleshooting — especially for config.vsh.yaml syntax,
    engine/runner setup, or Viash CLI usage.

    Args:
        query (str): Question about Viash, Docker, or Nextflow.

    Returns:
        str: Relevant content fetched live from viash.io.
    """
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

        scored = sorted(
            urls,
            key=lambda u: sum(1 for kw in keywords if kw in u.lower()),
            reverse=True,
        )

        session = requests.Session()
        session.headers["User-Agent"] = "Mozilla/5.0 (viash-docs-searcher)"
        parts = []
        for url in scored:
            text = _fetch_page(session, url)
            if text:
                parts.append(f"[Source: {url}]\n{text[:3000]}")
            if len(parts) >= 3:
                break

        result = "\n\n---\n\n".join(parts) if parts else "No relevant viash.io pages found."
        print(result)
        return result
    except Exception as e:
        err = f"Error searching viash.io: {e}"
        print(err)
        return err


# ─── Main initializer ─────────────────────────────────────────────────────────

def initialize_agent(
    model_id: str = None,
    llm_provider: str = None,
):
    """Return an A1 orchestrator with one tool and know-how docs injected."""

    load_dotenv()

    provider = (llm_provider or LLM_PROVIDER).lower()
    model = model_id or MODEL_ID
    source = _SOURCE_MAP.get(provider)
    if source is None:
        raise ValueError(f"Unsupported provider '{provider}'. Choose 'openai', 'anthropic', or 'openrouter'.")

    if provider == "openai" and not os.environ.get("OPENAI_API_KEY"):
        _missing_key_error("OPENAI_API_KEY")
    elif provider == "anthropic" and not os.environ.get("ANTHROPIC_API_KEY"):
        _missing_key_error("ANTHROPIC_API_KEY")
    elif provider == "openrouter":
        if not os.environ.get("OPENROUTER_API_KEY"):
            _missing_key_error("OPENROUTER_API_KEY")
        # Redirect the OpenAI SDK (used by Biomni A1) to OpenRouter's endpoint
        os.environ["OPENAI_API_KEY"] = os.environ["OPENROUTER_API_KEY"]
        os.environ["OPENAI_BASE_URL"] = "https://openrouter.ai/api/v1"

    agentic_root = Path(__file__).parent.parent

    print("Creating A1 orchestrator...")
    from biomni.agent.a1 import A1

    a1 = A1(
        path=str(agentic_root / ".biomni_data"),
        llm=model,
        source=source,
        use_tool_retriever=False,
        expected_data_lake_files=[],
        max_tokens=8192,
    )

    # Inject know-how docs into A1's know_how_loader so configure() includes them
    know_how_dir = agentic_root / "know_how"
    if know_how_dir.exists():
        n_loaded = 0
        for md_path in sorted(know_how_dir.glob("*.md")):
            try:
                content = md_path.read_text(encoding="utf-8")
                desc = ""
                for line in content.splitlines():
                    if line.startswith("**Short Description**"):
                        desc = line.split(":", 1)[-1].strip().strip("*").strip()
                        break
                a1.know_how_loader.add_custom_document(
                    doc_id=f"grn_{md_path.stem}",
                    name=md_path.stem.replace("_", " ").title(),
                    description=desc or md_path.stem,
                    content=content,
                )
                n_loaded += 1
            except Exception as e:
                print(f"Warning: could not load {md_path.name}: {e}")
        print(f"Injected {n_loaded} know-how doc(s) into A1")

    # Register search_viash by injecting schema directly — avoids an LLM call
    # that function_to_api_schema() would otherwise make during add_tool().
    print("Registering tools...")
    _search_viash_schema = {
        "name": "search_viash",
        "module": "agent",
        "description": (
            "Search the official Viash documentation (viash.io/guide and viash.io/reference) "
            "for information about Viash component configuration, Docker engine setup, "
            "Nextflow runner, argument types, resources, test benches, and CLI commands."
        ),
        "parameters": {
            "type": "object",
            "properties": {
                "query": {
                    "type": "string",
                    "description": "Question about Viash, Docker, or Nextflow.",
                }
            },
            "required": ["query"],
        },
        "required_parameters": ["query"],
    }
    import sys as _sys
    import types as _types
    _mod = _types.ModuleType("agent")
    _mod.search_viash = search_viash
    _sys.modules["agent"] = _mod
    if a1.use_tool_retriever:
        a1.tool_registry.register_tool(_search_viash_schema)
        print("Tool 'search_viash' registered via direct schema injection")
    else:
        # Inject directly into module2api — avoids the LLM call in add_tool()
        if "agent" not in a1.module2api:
            a1.module2api["agent"] = []
        a1.module2api["agent"].append(_search_viash_schema)
        print("Tool 'search_viash' injected directly into module2api")

    # Prepend routing instructions so A1 consults know-how context before tools.
    routing_prefix = (
        "## GRN Benchmark Assistant — Routing Rules\n\n"
        "You are an assistant for the `task_grn_inference` GRN benchmarking framework.\n"
        "Your system prompt already contains know-how documents covering: available GRN methods, "
        "datasets, evaluation metrics, output format, how to add a new method, and Viash commands.\n\n"
        "ROUTING RULES (follow strictly):\n"
        "1. For questions about GRN methods, datasets, evaluation, adding a method, or output format: "
        "answer DIRECTLY from the know-how context in your system prompt. Do NOT call search_viash.\n"
        "2. Only call search_viash when the user asks about Viash/Docker/Nextflow syntax or config "
        "details that are NOT already covered in the know-how docs.\n"
        "3. Never hallucinate. If the know-how docs don't cover something, say so.\n\n"
    )
    a1.system_prompt = routing_prefix + a1.system_prompt

    print("A1 agent ready!")
    return a1
