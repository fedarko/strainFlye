#! /usr/bin/env bash
# This is a short script that automatically overwrites the bottom lines of
# README.md with an up-to-date version of strainFlye's help text.
# It uses the unique word "STARTDOCS" in the README file as a sign for where
# to begin writing.

TMPFILENAME="_tempfile_pleasedontclobberme.txt"

# Maybe some year I'll remember how bash if statements work without having to
# look them up. 2021, apparently, is not that year.
# https://stackoverflow.com/a/638980
if [ -f $TMPFILENAME ]; then
    echo "$TMPFILENAME already exists? Don't wanna delete your stuff, so quitting."
fi

# Extract the line number the word STARTDOCS occurs at in the README:
# https://askubuntu.com/a/584410 and https://askubuntu.com/questions/77156/how-to-get-line-number-from-grep#comment1734714_584410
DOCSLINENUM=`grep -n "STARTDOCS" < README.md | cut -f 1 -d :`

# Store the stuff in README.md we want to keep in a temporary file.
# Unfortunately, if we try to redirect straight back to README.md, we end up
# with an empty file: see https://stackoverflow.com/a/339941. Hence why we use
# a temporary file.
head -n $DOCSLINENUM README.md > $TMPFILENAME

# Now, clobber README.md so that it only includes the lines up to and including
# the STARTDOCS line.
cat $TMPFILENAME > README.md

# Now we can safely add in the help text (and surrounding ``` code block
# signs).
echo '```' >> README.md
strainFlye -h >> README.md
echo '```' >> README.md

# Remove the temporary file. Note that we have already checked at the start of
# this file to make sure that this file didn't already exist, so there
# shouldn't be a danger of accidentally removing a real file that happens to
# exist with this name in this directory. (... Unless someone would go to the
# trouble of trying to exploit a race condition and add text to this file after
# the above check but before this line. I don't think we need to worry about
# that here.)
rm $TMPFILENAME
