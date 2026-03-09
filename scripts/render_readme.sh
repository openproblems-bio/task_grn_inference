set -e

common/scripts/create_task_readme --input src/api

# Inject info.readme from _viash.yaml into the generated README.md,
# placing it between the "Repository:" line and "## Description".
python3 - <<'EOF'
import re

# Parse the literal block scalar under info.readme without third-party deps
with open("_viash.yaml") as f:
    lines = f.readlines()

in_readme = False
readme_indent = None
readme_lines = []

for line in lines:
    if re.match(r'\s+readme:\s*\|\s*$', line):
        in_readme = True
        readme_indent = None
        continue
    if in_readme:
        if line.strip() == "":
            readme_lines.append("")
            continue
        current_indent = len(line) - len(line.lstrip())
        if readme_indent is None:
            readme_indent = current_indent
        if current_indent < readme_indent:
            break
        readme_lines.append(line[readme_indent:].rstrip())

readme_insert = "\n".join(readme_lines).strip()
if not readme_insert:
    raise SystemExit("No info.readme found in _viash.yaml")

with open("README.md") as f:
    content = f.read()

# Insert after the "Repository: ..." line and before "## Description"
pattern = r'(Repository:.*?\n)(\s*\n## Description)'
replacement = r'\1\n' + readme_insert + r'\n\2'
new_content, n = re.subn(pattern, replacement, content, flags=re.DOTALL)

if n == 0:
    raise SystemExit("Could not find insertion point in README.md")

with open("README.md", "w") as f:
    f.write(new_content)

print("Injected info.readme into README.md")
EOF