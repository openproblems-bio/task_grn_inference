# GeneRNBI Agentic Assistant

A conversational agent that helps users navigate and extend the GeneRNBI benchmarking framework. It uses a [LlamaIndex ReAct](https://docs.llamaindex.ai/en/stable/examples/agent/react_agent/) reasoning loop and answers questions by searching the GeneRNBI documentation.

## Structure

```
agentic/
├── src/
│   ├── run.py           # Entry point (interactive or single-query)
│   ├── agent.py         # ReAct agent definition
│   ├── rag_builder.py   # Builds/loads VectorStoreIndex over docs/
│   ├── config.py        # LLM provider and model selection
│   └── instructions.txt # System prompt
├── docs/                # Documentation files indexed by the agent
│   ├── README.md
│   └── METRIC_EVALUATION_README.md
├── run.sh               # Convenience shell script
└── .env.template        # API key template (do NOT commit a filled .env)
```

## Setup

### 1. Conda environment

The agent requires the `llamaindex` conda environment:

```bash
conda activate llamaindex
```

If you don't have it, install dependencies into any Python 3.10+ environment:

```bash
pip install llama-index-core llama-index-llms-openai llama-index-llms-anthropic python-dotenv
```

### 2. API key

> **Keys must never be committed to the repo.** The `.env` file is gitignored.

First, get your API key:

| Provider | Where to get it |
|----------|----------------|
| **OpenAI** | [platform.openai.com/api-keys](https://platform.openai.com/api-keys) |
| **Anthropic** | [console.anthropic.com/settings/keys](https://console.anthropic.com/settings/keys) |

Then provide it via **one** of these methods:

**Option A — inline, one-off run:**
```bash
OPENAI_API_KEY=sk-... bash run.sh
# Anthropic:
ANTHROPIC_API_KEY=sk-ant-... bash run.sh
```

**Option B — export for the current session:**
```bash
export OPENAI_API_KEY=sk-...
bash run.sh
```

**Option C — permanent (recommended): add to `~/.bashrc` or `~/.zshrc`:**
```bash
echo 'export OPENAI_API_KEY=sk-...' >> ~/.bashrc
source ~/.bashrc
```

**Option D — `.env` file (local only):**
```bash
cp .env.template .env
# open .env and fill in your key
```

If the key is missing, the agent will print a clear error with these instructions.

### 3. Select LLM provider / model

Edit `src/config.py`:
```python
LLM_PROVIDER = "openai"    # or "anthropic"
MODEL_ID     = "gpt-4o"    # or "claude-sonnet-4-6", etc.
```

## Running

**Interactive mode:**
```bash
bash run.sh
# or directly:
cd src && python run.py
```

**Single-query mode:**
```bash
bash run.sh "What layer name does the evaluation script expect?"
# or:
cd src && python run.py "How do I integrate a new GRN method?"
```

## Extending the knowledge base

Drop any `.md`, `.txt`, or `.pdf` files into `docs/` and delete `.index_cache/` to trigger a rebuild on next run:

```bash
rm -rf .index_cache/
bash run.sh
```
