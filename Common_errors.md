# Common Errors
### 1.1 Perl: Can't locate NAME.pm

If you try to run a perl script and get an error message that starts like this:

```
Can't locate List/MoreUtils.pm in @INC

```

This means you are missing a required perl module to run the script. Try the following:

```
ml perl

cpan [package to install]
```

For example, for the above error you would run:

```
cpan List::MoreUtils
```


### 1.2 Sherlock modules: These module(s) exist but cannot be loaded

This error after trying "module load" means that the module cannot be loaded until another module is loaded first.

In the case of commonly used programs in our lab usually this mean you need to load the biology modules first,
i.e.:

```
module load biology
```


