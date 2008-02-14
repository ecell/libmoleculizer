
make clean

# This has to happen before the first dose/response suite is run by the
# optimizer, so that it will exist when the optimizer wants to journal optimal
# values into it.

# moleculizer-variations is removed by preen, but not by clean.
make moleculizer-variations

# This does the actual dose/response suite to evaluate transinformation.
make
