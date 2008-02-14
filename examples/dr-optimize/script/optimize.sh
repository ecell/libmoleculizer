

# Note that if the min-size option is not given, then no optimization
# is done, and only the "empirical" transinformation is used instead.

procOpt \
    --proc-name ./script/remake.sh \
    --proc-input-file-name ./substitution-values \
    --proc-output-file-name ./dose-response.out/tddrs/active-kinases--active-kin3--transinfo \
    --step-size-file-name ./substitution-steps \
    --journal-file-name ./optimization-journal \
    --journal-proc-name ./script/journal-optimize-cycle.pl \
    --min-size 0.001

