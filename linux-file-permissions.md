# Linux File Permissions

> Based on the Red Hat blog post [Linux file permissions explained](https://www.redhat.com/en/blog/linux-file-permissions-explained).

## 1. Why does it matter?

Set your file permissions so that other members of the lab group can read your files. This matters because:

1. Other lab members might need access to replicate your code.
2. Lab member folders are backed up from **Oak** to **Elm**, and the backup fails if some of your files are not accessible.
3. It is especially **important** when you are about to leave the lab group.
4. It matters most for files in the **Oak** folder. Files in your own home or scratch directory don't need to be shared.

## 2. Who permissions apply to

Each file and directory has three sets of permissions:

| Class | Symbol | Meaning |
|-------|--------|---------|
| User (owner) | `u` | The account that owns the file (you) |
| Group | `g` | Members of the file's group (`schumerlab`) |
| Others | `o` | Everyone else |
| All | `a` | All three of the above (used with `chmod`) |

## 3. Permission types

In general, everyone in the lab group should be able to **read** (and, where needed, **execute**) your files. Only you, the owner, should have **write** permission, so no one can accidentally change your scripts.

| Permission | Symbol | On a file | On a directory |
|------------|--------|-----------|----------------|
| Read | `r` | View the contents | List the entries (`ls`) |
| Write | `w` | Change the contents | Create, delete, or rename entries |
| Execute | `x` | Run it as a program/script | Enter it (`cd`) and reach files inside |

> **Note:** Directories need `x` as well as `r`. Without `x`, lab members can't open the folder or reach the files inside, even if those files are readable.

## 4. Checking permissions with `ls -l`

To see the current permissions of a file:

```bash
$ ls -l script.sh
-rwxr-xr-- 1 weiluwl oak_schumer 512 Nov 10 09:00 script.sh
```

Breaking down `-rwxr-xr--`:

```
-    rwx    r-x    r--
│     │      │      └── others: read
│     │      └───────── group:  read, execute
│     └──────────────── user:   read, write, execute
└────────────────────── file type
```

The owner (`weiluwl`) and group (`oak_schumer`) are listed right after the permissions. Make sure the group is `oak_schumer`; if it isn't, group permissions won't help lab members.

Common file type characters: `-` regular file, `d` directory, `l` symbolic link.

## 5. Changing permissions with `chmod`

Two common cases:

1. **Making a script executable.** New text files are usually read/write only. To run a shell script you created, add execute permission.
2. **Sharing files with the lab.** Make sure the group can read your files and open your directories.

### Symbolic mode

Format: `chmod [who][+|-|=][permissions] file`

```bash
chmod u+x script.sh    # let yourself run the script
chmod g+r file.txt     # let the lab group read the file
chmod go-w file.txt    # remove write access from group and others
```

### Octal (numeric) mode

Numeric mode is often easier than typing out each permission. Each permission has a value; add them up for each class (user, group, others):

| Permission | Value |
|------------|-------|
| `r` | 4 |
| `w` | 2 |
| `x` | 1 |
| none | 0 |

| Octal | Symbolic | Meaning |
|-------|----------|---------|
| 7 | `rwx` | read + write + execute |
| 6 | `rw-` | read + write |
| 5 | `r-x` | read + execute |
| 4 | `r--` | read only |
| 0 | `---` | no access |

```bash
chmod 755 script.sh   # rwxr-xr-x — scripts and directories (recommended)
chmod 644 file.txt    # rw-r--r-- — regular files (recommended)

chmod 600 id_rsa      # rw------- — private to owner
chmod 700 ~/private   # rwx------ — owner-only directory
```

**Recommended:** `755` for directories and scripts, `644` for regular files.

### Recursive

By default, `chmod` changes only the file or folder you name, not what's inside it. Use `-R` to apply changes to everything inside a directory.

Avoid `chmod -R 755` on a whole folder: it makes every file executable. Use this instead:

```bash
chmod -R g+rX /path/to/dir
```

Capital `X` adds execute only to directories (and files that are already executable). This gives the lab read access to everything and lets them open every subfolder.

## 6. Quick cheat sheet

| Task | Command |
|------|---------|
| View permissions | `ls -l file` |
| Make a script executable | `chmod +x script.sh` |
| Standard file | `chmod 644 file` |
| Standard directory / script | `chmod 755 dir` |
| Share a whole folder with the lab | `chmod -R go+rX dir` |
| Fix the group of a folder | `chgrp -R schumerlab dir` |
