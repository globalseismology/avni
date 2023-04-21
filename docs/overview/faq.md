(faq)=
```{eval-rst}
.. include:: ../links.inc
.. include:: ../buttons.inc
```

# Frequently Asked Questions (FAQ)

```{eval-rst}
.. highlight:: python
```
## General queries

### Is it free?

Yes! The master branch of our repository that hosts the client-side Python codes is synchronized for open-source access from various locations on the internet - pypi, conda-forge - and is hosted in a public Github account. This repository can be used for free under the terms of [our license](../getting-started/cite).


### Is all your code open-source?

Our public Github pure-Python repository is provided open-source with the GNU GPL v3 [license](../getting-started/cite). In order to encourage involvement by the relevant domain experts without the necessary overhead of immediately catering to requests in an open-source environment, we keep a portion of the hard-to-compile Fortran, C and other Python routines private to the AVNI development team. Feel free to raise an [issue](https://github.com/globalseismology/avni/issues) that you can help with or write to use at **avni@globalseismology.org** if you want to get involved.

### How can we get involved?

Become a tester or contributor! Please try out our codes in various applications and let us know. Fork our public repository, contribute code and raise [issues or requests](https://github.com/globalseismology/avni/issues). If you want to be a co-developer, please request access through our [Website](http://avni.globalseismology.org).

### I am a domain expert in geosciences. Why should I get involved?

Because you will have fun and work with other domain experts! We have the necessary infrastructure in place for benchmarking, developing and testing algorithms. Your work will be preserved for other colleagues to build on and benchmark against. You will leverage legacy codes spanning multiple decades, modified and optimized for modern infrastructure. You will work in a copyrighted, private branch or repository that will not be released publicly or used in other projects or branches without your written approval. Given other scientific commitments, you will avoid the necessary overhead of immediately catering to requests inherent in other standard open-source environments.

### How do I cite AVNI?

See {ref}`cite`.

### I'm not sure how to do *X* analysis step with my *Y* data...

Knowing "the right thing" to do with geoscience data is challenging. We use
the [AVNI Forum](https://github.com/globalseismology/avni/discussions) to discuss analysis strategies for different kinds of
data. It's worth searching the archives to see if there have been relevant
discussions in the past, but don't hesitate to ask a new question if the answer
isn't out there already.

## Installation

### Help! I can't get Python and AVNI working!

Check out our {ref}`installation instructions <installers>`.

### I still can't get it to work!

See {ref}`help`.

### What are the system requirements?

AVNI has been tested on the following platforms: Linux (Ubuntu, Redhat, CentOS), MacOS and Windows. It requires a Python installation with versions 3.3 or above.


### Python runs on macOS extremely slow even on simple commands!

Python uses some backends that interfere with the macOS energy saver when
using an IDE such as Spyder or PyCharm. To test it, import `time` and run:

```
start = time.time(); time.sleep(0.0005); print(time.time() - start)
```

If it takes several seconds you can either:

- Install the module `appnope` and run in your script:

  ```
  import appnope
  appnope.nope()
  ```

- Change the configuration defaults by running in your terminal:

  ```console
  $ defaults write org.python.python NSAppSleepDisabled -bool YES
  ```

### I think I found a bug, what do I do?

When you encounter an error message or unexpected results, it can be hard to
tell whether it happened because of a bug in AVNI, a mistake in user
code, a corrupted data file, or irregularities in the data itself. Your first
step when asking for help should be the
[AVNI Forum](https://github.com/globalseismology/avni/discussions), not GitHub. This bears
repeating: *the GitHub issue tracker is not for usage help* — it is for
software bugs, feature requests, and improvements to documentation. If you
open an issue that contains only a usage question, we will close the issue and
direct you to the forum. If you're pretty sure the problem you've encountered
is a software bug (not bad data or user error):

- Make sure you're using [the most current version]. You can check it locally
  at a shell prompt with:

  ```console
  $ avni sys_info
  ```

  which will also give you version info about important AVNI
  dependencies.

- If you're already on the most current version, if possible try using
  {ref}`the latest development version <installing_main>`, as the bug may
  have been fixed already since the latest release. If you can't try the latest
  development version, search the GitHub issues page to see if the problem has
  already been reported and/or fixed.

- Please provide a
  link to a small, anonymized portion of your data that does yield the error.

If the problem persists, [open a new issue](https://github.com/globalseismology/avni/issues/new)
and include the *smallest possible* code sample that replicates the error
you're seeing. Paste the code sample into the issue, with a line containing
three backticks (\`\`\`) above and below the lines of code. This
minimal working example should be self-contained, which means that
AVNI contributors should be able to copy and paste the provided snippet
and replicate the bug on their own computers.


## API and Privacy

### Why do you provide APIs?

We provide application programming interfaces (APIs) that interface with heavy, legacy codes hosted on our servers so that AVNI installation remains light to serve various applications. We care deeply about facilitating science by reducing the time a typical researcher or student spends installing complex dependencies and debugging code. Some users may not have the necessary infrastructure to deal with big data and our scientific codes. It is very difficult for us to spend time testing our server-side codes across platforms and transfer terabytes of data to a client computer.

## Is there a cap on the number of calls to the API?

The number of calls are capped at 1000 per day per user. This is limited by the current hardware and as the project evolves, we hope to increase the limit. If you find AVNI useful or want more features, please let the funding agencies know or leave a public [comment](https://github.com/geodynamics/avni/issues).