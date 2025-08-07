bash_file="scripts/prior/run_process_data.sh"

jid=$(sbatch $bash_file | awk '{print $4}'); \
out_file="logs/${jid}.out"; \
err_file="logs/${jid}.err"; \
while [ ! -f "$out_file" ] || [ ! -f "$err_file" ]; do sleep 1; done; \
echo "$out_file" is running