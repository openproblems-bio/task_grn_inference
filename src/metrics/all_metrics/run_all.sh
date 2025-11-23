run_type='sbatch'

metrics=(
    rc_tf_act 
    regression
    tf_recovery
    ar
    rc_tf_act
    tf_binding
    vc
    ws_distance
    sem
    gs_recovery
)
for metric in "${metrics[@]}"; do
    echo "Running metric: ${metric}"
    ${run_type} src/metrics/${metric}/run_local.sh 
    echo "----------------------------------------"
done

