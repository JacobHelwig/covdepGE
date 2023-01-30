cd ../../../

# save dir
save_dir="./experiments/discrete/cdge"

# experiment name
exp=disc_cov_dep

# sample size and signal strength


# number of trials
ntrials=50

# data parameters
n=100
dims=( 11 31 51 )
n1s=( 50 80 )
lambdas=( 3 15 )

# models to skip
skips="c('mgm','varbvs')"

for p in ${dims[@]}; do
  for n1 in ${n1s[@]}; do
    for lambda in ${lambdas[@]}; do
      n2=$(( n - n1 ))
      now=$(date '+%Y-%m-%d_%H:%M')
      out_name=$save_dir/$exp\_ntrials$ntrials\_p$p\_n1\_$n1\_n2\_$n2\_lambda$lambda\_$now.Rout
      printf "\n\nRunning $out_name...\n"
      SECONDS=0
      R CMD BATCH --no-save --no-restore "--args \
      save_dir='$save_dir' \
      experiment='$exp' \
      p=$p \
      n1=$n1 \
      n2=$n2 \
      lambda=$lambda \
      n_trials=$ntrials \
      skips=$skips" \
      main.R $out_name
      duration=$SECONDS
      echo "Done. $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed"
    done
  done
done
