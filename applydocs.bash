doxygen
git config merge.renameLimit 999999
git add -f docs
git stash 
git checkout gh-pages
git stash apply
grep -lr '<<<<<<< ' . | xargs git checkout --theirs
git add -u .
git config --unset merge.renameLimit
