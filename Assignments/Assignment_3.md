# Assignment 3 [*60 marks*]

## Note that failure to include a coversheet will incur a 5 mark penalty.

## Assignment overview [*8 marks for general format*]

This assignment is based on the Local version control and GitHub repository practicals. 

It is highly recommended that you attend "GitHub Advanced" Practical as we will be walking through much of this assignment (generated to have different answers)

The repository at https://github.com/davmlaw/2026_assignment_3_version_control is a simulation of an open source project with bugs.

This assignment involves performing actions on GitHub [*8 marks total for general format*]:

* Issue for XXXXXXX bug [*2 marks*]
* Issue for "UPPERCASE" bug [*2 marks*]
* Fork repo [*2 marks*]
* Pull request for fixing the issues [*2 marks*] 

These have timestamps on them, and these will be checked against the submission date.

Answers to theory questions should be submitted to MyUni as a file (plain text or PDF)

## Practical questions [*42 marks*] 

Please read the instructions in the initial README of the repository, and only modify the files assigned to you.

There are two kinds of changes/bugs:

* insertion of XXXXXXX
* changing some lower case words to UPPERCASE

These should be treated as two separate incidents, with different severities, and will be handled separately.

1. Raise an issue in the [repo](https://github.com/davmlaw/2026_assignment_3_version_control) for the XXXXXXX bug
    + Giving a good title succinctly stating the problem and identifying the file (such that it is clear why it is distinct from other issues) [*2 points*]
    + Thank the maintainer and use a polite, kind and non-accusatory tone  [*2 points*]
    + Provide the commit message of the commit that introduced the error to your file [*2 points*]
    + Write an expected / actual result for the lines affected, only mentioning the high severity bug (XXXXXXX) [*4 points*]
    + Add a comment that you plan on submitting a pull request [*2 points*]

(total: 12 points)

2. Raise an issue in the [repo](https://github.com/davmlaw/2026_assignment_3_version_control) for the UPPERCASE bug
    + Giving a good title succinctly stating the problem and identifying the file (such that it is clear why it is distinct from other issues) [*2 points*]
    + Provide the date and Git hash of the commit that introduced the error to your file [*2 points*]
    + Write an expected / actual result for the lines affected, only mentioning the low severity bug (UPPERCASE) [*4 points*]

(total: 8 points)

3. Fork the repository to your personal account.
    + *Make sure you untick* the "Copy the main branch only" - you want the other "stable" branch to also be copied

4. Fix the XXXXXXX issue, in the "main" and "stable" branches
    + You will need to clone your private fork of the repo to your local machine.
    + Edit your data/ file, and remove the XXXXXXX ONLY (ie DO NOT change the UPPERCASE at the same time). There should only be 1 space between the words after removing the XXXXXXX [*2 points*]
    + Edit your docs/ file, which is in the format of a [changelog](https://keepachangelog.com/en/1.0.0/) adding a CHANGED entry for your fix, that links to the high severity issue above [*2 points*]
    + Write a descriptive commit message, which references the high severity issue above [*2 points*]
    + This bug needs to be fixed in the "stable" branch as well, so switch to that branch and "cherry-pick" the commit you just made [*6 points*]  
    + Push both the "main" and "stable" branches to your repo

(total: 12 points)

5. Fix the UPPERCASE issue, in the "main" branch ONLY
    + Edit your data/ file, and remove the low severity issue. [*2 points*]
    + Edit your docs/ file, which is in the format of a [changelog](https://keepachangelog.com/en/1.0.0/) adding a CHANGED entry for your fix, that links to the low severity issue above [*2 points*]
    + Write a descriptive commit message, which references the low severity issue above [*2 points*]
    + Push the main branch to your repo

(total: 6 points)

6. Submit a pull request from your main branch. This pull request should be for both issues/fixes. In real life we'd do a separate pull request for each fix, but for this assignment we'll do them both in 1 PR.

    + Comments should link to both issues and give a brief description of the change [*4 points*] 

(total: 4 points)
    
## Theory questions [*10 marks*]

1. How do you tell Git to stop tracking / ignore some files. Name 2 types of files you might want to ignore and why. [*2 points*]

2. What is the staging area in Git, and why is it used? [*2 points*]

3. Why would you want a long lasting branch (eg "stable" above). [*2 points*]
 
4. Why would you want a short lived one (ie a "dev" or "feature" branch) [*2 points*]

5. Why is it called a pull request? [*2 points*]

### Resources

- [How to write a good bug report](https://musescore.org/en/node/309537)
- [Issues - Actual vs Expected](https://medium.com/we-are-testers/chapter-2-how-to-write-useful-actual-and-expected-results-details-in-your-bug-report-10b83e5aaa75)
- [Git - branches in a nutshell](https://git-scm.com/book/en/v2/Git-Branching-Branches-in-a-Nutshell)
- [GitHub - about branches](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-branches)
- [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
- [GitHub - creating a pull request from a fork](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork)
