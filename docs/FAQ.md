Frequently Asked Questions
--------------------------

  1. [General](#general)  
    * [Is it free?](#is-it-free)  
    * [Is all your code open-source?](#is-all-your-code-open-source)  
    * [How can we get involved?](#how-can-we-get-involved)  
    * [I am a domain expert in geosciences. Why should I get involved?](#i-am-a-domain-expert-in-geosciences-why-should-i-get-involved)  
  2. [Installation](#installation)  
    * [What are the system requirements for REM3D?](#what-are-the-system-requirements-for-rem3d)  
  3. [Privacy](#privacy)  
    * [Why do you provide APIs?](#why_do_you_provide_apis)  
    * [Is there a cap on the number of calls to the API?](#is-there-a-cap-on-the-number-of-calls-to-the-api) 


General
-------

Is it free?
-----------

Yes! The master branch of our repository that hosts the client-side Python codes is synchronized for open-source access from various locations on the internet - pypi, conda-forge - and is hosted in a public Github account. This repository can be used for free under the GNU GPL v3 [LICENSE](../LICENSE).

Is all your code open-source?
----------------------------

Our public Github pure-Python repository is provided open-source with the GNU GPL v3 [LICENSE](../LICENSE). In order to encourage involvement by the relevant domain experts without the necessary overhead of immediately catering to requests in an open-source environment, we keep a portion of the hard-to-compile Fortran, C and other Python routines private to the REM3D development team. Feel free to raise an issue [issue](https://github.com/globalseismology/rem3d/issues) that you can help with or write to use at **info@rem3d.org** if you want to get involved.

How can we get involved?
------------------------

Become a tester or contributor! Please try out our codes in various applications and let us know. Fork our public repository, contribute code and raise [issues or requests](https://github.com/globalseismology/rem3d/issues). If you want to be a co-developer, please request access through our [Website](http://rem3d.org/join-us/github).

I am a domain expert in geosciences. Why should I get involved?
---------------------------------------------------------------

Because you will have fun and work with other domain experts! We have the necessary infrastructure in place for benchmarking, developing and testing algorithms. Your work will be preserved for other colleagues to build on and benchmark against. You will leverage legacy codes spanning multiple decades, modified and optimized for modern infrastructure. You will work in a copyrighted, private branch or repository that will not be released publicly or used in other projects or branches without your written approval. Given other scientific commitments, you will avoid the necessary overhead of immediately catering to requests inherent in other standard open-source environments.

Installation
------------

What are the system requirements?
--------------------------------

REM3D has been tested on the following platforms: Linux (Ubuntu, Redhat, CentOS), MacOS and Windows. It requires a python installation with versions 3.3 and above.

Privacy
-------

Why do you provide APIs?
------------------------

We provide application programming interfaces (APIs) that interface with heavy, legacy codes hosted on our servers so that REM3D installation remains light to serve various applications. We care deeply about facilitating science by reducing the time a typical researcher or student spends installing complex dependencies and debugging code. Some users may not have the necessary infrastructure to deal with big data and our scientific codes. It is very difficult for us to spend time testing our server-side codes across platforms and transfer terabytes of data to a client computer. 

Is there a cap on the number of calls to the API?
-------------------------------------------------

The number of calls are capped at 5000 per day per public user. This is limited by the current hardware and as the project evolves, we hope to increase the limit. If you find REM3D useful or want more features, please let the funding agencies know or leave a public [comment](https://github.com/globalseismology/rem3d/issues).