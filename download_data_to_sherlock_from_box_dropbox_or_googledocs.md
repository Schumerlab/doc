# Download data to Sherlock from box, dropbox or googledocs

## Contents
- 1. rclone can download data to Sherlock
    - 1.1 [configure download from dropbox using rclone](#configure-download-from-dropbox-using-rclone)
    - 1.2 [configure download from Box using rclone](#configure-download-from-box-using-rclone)
    - 1.3 [Other access options with rclone](#other-access-options-with-rclone)

## rclone can download data to Sherlock
rclone is super useful and convenient and especially for large files, but is a little annoying for first time use.

### configure download from dropbox using rclone
1) Make sure you are logged into dropbox in your browser
2) Follow this link (https://web.archive.org/web/20230602223208/https://rclone.org/downloads/) to download the appropriate desktop version of rclone(usually 32 bit macOS or 32 bit Windows):
    - unzip the archive
    - open a terminal window on your computer
    - navigate to this folder in your Downloads. For example, on a macOS:

    `cd Downloads/rclone-v1.46-osx-386`

3) Now, open another terminal window, leaving this one open, and login to sherlock and type:

    `ml load system rclone`

4) Configure rclone:

    `rclone config`

    
You will be prompted:

```
No remotes found - make a new one
n) New remote
s) Set configuration password
q) Quit config
n/s/q>
```

Respond with: `n`

You will be prompted to provide a name, you can provide any name - e.g.:

`name> remote_dropbox`

You will be prompted to select what kind of storage you are configuring to, enter 7 (with rclone v1.62.2 64bit macOS, enter 13):

`Storage>7`

You will be prompted to select a client_id, leave this blank:

`client_id>`

You will be prompted to select a client_secret, leave this blank:

`client_secret>`

You will be asked if you want to edit the advanced config, enter n:

`Edit advanced config? (y/n)>n`

You will be asked if you want to auto config, say no!:

```
Use auto config?
* Say Y if not sure
* Say N if you are working on a remote or headless machine
```

5) Now, switch to your desktop terminal in the rclone directory. Past the following there:

`./rclone authorize "dropbox"`

This will open a browser window and ask you to authorize, after you do so, a code will appear in your desktop
terminal. Copy between:

Paste the following into your remote machine *---> code will be here <---* End paste

6) Enter the copied text into the Sherlock terminal. Then respond the following prompt with y:

```
y) Yes this is OK
e) Edit this remote
d) Delete this remote
y/e/d> y
```

7) This will complete the configuration but it's recommended that you add a password:
remote_dropbox dropbox

```
e) Edit existing remote
n) New remote
d) Delete remote
r) Rename remote
c) Copy remote
s) Set configuration password
q) Quit config
e/n/d/r/c/s/q> s
```

When prompted respond with a to add a password:
```

a) Add Password
q) Quit to main menu
a/q> a
```

follow the prompts to add your password

8) You can then view your available files and folders in dropbox:

`rclone lsd remote_dropbox:`

or copy or sync them:


`rclone copy remote_dropbox:Schumer_lab_resources/README_depositing_data.txt
./`

for full list of options, try:

`rclone --help`

### configure download from Box using rclone
repeat these steps for Box, except give it a different name (e.g. remote_box), and select 5 instead of 7 (with rclone v1.62.2 64bit macOS, enter 8)

access with:

`rclone lsd remote_box:`

### Other access options with rclone
See here for options of storage endpoints compatible with rclone (https://web.archive.org/web/20230602223208/https://rclone.org/)