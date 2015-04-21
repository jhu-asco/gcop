doxygen
git add -f docs
git stash 
git checkout gh-pages
git stash apply
grep -lr '<<<<<<<' . | xargs git checkout --theirs
git add -u .
