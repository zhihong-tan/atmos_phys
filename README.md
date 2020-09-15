# atmos_phys
This is the repository that combines the atmosphere physics and chemistry
## Instructions on merging in code from atmos_param or atmos_shared
In order to commit changes that currently exist on atmos_param or atmos_shared,
you first have to commit and push the changes to a branch in the repository you are 
working in. Next you will clone atmos_phys.  After that, you will add the original 
repository to as a remote in your atmos_phys directory.  You can then fetch from the
original repository and merge in the branch you were working on.

Below is an example set of instructions to accomplish this.  This assumes that you 
cloned atmos_param.
```bash
git clone https://gitlab.gfdl.noaa.gov/fms/atmos_phys.git #clone atmos_param
cd atmos_phys 
pushd ../atmos_param/ #Set up a directory stack with the two repo directories
git checkout user/ter/exampleBranch #check out your user branch
git fetch #fetch any changes that may have occurred in atmos_param since you branched
git pull origin/master #merge in the master (update your code base)
## The following 3 lines can be skipped if you already have committed changes ##
touch test.file #Make your changes (DO NOT COPY THIS, it is just an example)
git add test.file #stage your changes
git commit -m "Save a test commit to user branch" #commit your changes
## The folowing can be skipped if you have already pushed your changes ##
git push -u origin user/ter/exampleBranch
pushd #change directory to the atmos_phys folder
git checkout -b user/ter/exampleBranch #create your branch in atmos_phys
git remote add param https://gitlab.gfdl.noaa.gov/fms/atmos_param.git #add atmos_param as a remote
git fetch param #fetch the contents of the atmos_param
git merge param/user/ter/exampleBranch #merge your branch from atmos_param into your branch on atmos_phys
git push -u origin user/ter/exampleBranch
```
A note on git remotes: 

The default name for a remote is **origin**.  When you clone a repository, the remote 
is automatically given the name **origin**.  When you add a remote, you can name it 
anything except origin.  In the above example, I named the remote **param**.  That's 
why I fetch from **param**.  To see what remotes you have set up, you can issue the 
command
```bash
git remote -v
```
