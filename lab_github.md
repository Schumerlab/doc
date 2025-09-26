## Lab github

Our lab GitHub is maintained here:
https://github.com/Schumerlab

On sherlock, our lab repositories are linked to our shared bin:

`/home/groups/schumer/shared_bin/Lab_shared_scripts`

### Tutorials in using github
Once you have a GitHub account and have done some tutorials, ask to get added to our lab GitHub
Tutorials:

https://product.hubspot.com/blog/git-and-github-tutorial-for-beginners
https://kbroman.org/github_tutorial/

### Important!
When adding scripts to lab repositories that other people may be working on (like Lab_shared_scripts), it is
important to avoid conflicts: Make sure to pull the latest changes before you push to avoid conflicts!!!

For example:

```
git add myscript.pl
git commit -m "my new script"
git pull
git push origin master
```

If you are working on a commonly used script that other people might be working on check the latest commit
before pulling to make sure you will not face merge conflicts.