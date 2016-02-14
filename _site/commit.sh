git add --all
echo -n "Enter commit description >" 
read text
git commit -m "$text"
git push -u origin master
