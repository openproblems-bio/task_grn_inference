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
├── singularity/
│   ├── agentic.def      # Singularity definition file
│   └── build_and_push.sh# Build image and push to S3
├── run.sh               # Convenience shell script
└── .env.template        # API key template (do NOT commit a filled .env)
```

## Setup
create conda env genernbi_agentic --->
#TODO: this needs fixing based on the new biomni requitements.


```bash
pip install llama-index-core llama-index-llms-openai llama-index-llms-anthropic python-dotenv
```

### 2. API key

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

## Container (Singularity / Apptainer)

A pre-built Singularity image is available on S3 and includes all required
dependencies (Biomni A1, LlamaIndex, OpenAI / Anthropic clients, HuggingFace
embeddings).

### Download

```bash
aws s3 cp s3://openproblems-data/resources/grn/images/agentic.sif . --no-sign-request
```

### Run

**Interactive mode:**
```bash
export OPENAI_API_KEY=sk-...       # or ANTHROPIC_API_KEY for Anthropic
singularity run \
    --bind /path/to/task_grn_inference/agentic:/agentic \
    agentic.sif
```

**Single-query mode:**
```bash
export OPENAI_API_KEY=sk-...
singularity run \
    --bind /path/to/task_grn_inference/agentic:/agentic \
    agentic.sif "How do I integrate a new GRN method?"
```

**Via `run.sh` (convenience wrapper):**
```bash
export OPENAI_API_KEY=sk-...
USE_SINGULARITY=1 bash run.sh "How do I integrate a new GRN method?"
# Override image path if needed:
SINGULARITY_IMAGE=/path/to/agentic.sif USE_SINGULARITY=1 bash run.sh
```

The `--bind` mount maps the `agentic/` directory to `/agentic` inside the
container so the agent can find its `docs/`, `data/`, and source files.

### Build the image yourself

```bash
cd agentic/singularity
bash build_and_push.sh              # build + upload to S3
bash build_and_push.sh --build-only # build only
```

---

## Extending the knowledge base

Drop any `.md`, `.txt`, or `.pdf` files into `docs/` and delete `.index_cache/` to trigger a rebuild on next run:

```bash
rm -rf .index_cache/
bash run.sh
```
