
for TARGET in $*
do
  make -q $TARGET > /dev/null 2>&1
  if test $? -eq "1"
  then
      echo $TARGET
  fi
done
	