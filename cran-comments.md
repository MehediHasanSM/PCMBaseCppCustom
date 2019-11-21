---
title: "NEWS about the PCMBaseCpp R-package"
author: "Venelin Mitov"
date: "21 November, 2019"
output: html_document
---

Disabled OPENMP compilation in an attempt to fix the failures on "Fedora clang" and "Solaris". Note that these errors could not be reproduced using devtools::check_rhub(). Following is the e-mail exchange with Prof. Brian Ripley:


Am Do., 14. Nov. 2019 um 17:36 Uhr schrieb Venelin Mitov <vmitov@gmail.com>:

Dear Prof. Ripley, 
Thanks for your message and for running the installation in debug mode. I've tested the package using devtools::check_rhub() and I only get notes for fedora-clang (screenshot below). 
 

From the trace you've sent me, the error with fedora-clang seems to occur at line 1853 of SPLITT.h, which reads
this->ptr_ = std::current_exception(); 
Is it possible that the standard C++ function std::current_excpetion() does call terminateon f edora-clang? I don't have any way to reproduce the issue or test fixes on this platform. As for the solaris, there seems to be an additional error in the computation, because the logged values are taking random signs, while they should always be non-negative. Again, I don't have means to reproduce this issue:
  what(): QuadraticPoly.h:InitNode:: The matrix V for node 15 (branch length=13.9031) is nearly singular or not positive definite; near 0 or negative eigenvalue found: -113.107; V.slice(i)(ki,ki)


Best regards, Dr. Venelin Mitov


Am Do., 14. Nov. 2019 um 12:37 Uhr schrieb Prof Brian Ripley <ripley@stats.ox.ac.uk>:
See https://cran.r-project.org/web/checks/check_results_PCMBaseCpp.html .

This version consistently fails on fedora-clang and Solaris.  Note that 
a C++ program calling 'terminate' is contrary to the CRAN policies.  And 
that is what is happening too with clang: running under gdb shows

Thread 1 "R" received signal SIGSEGV, Segmentation fault.
0x00007ffff56752c8 in typeinfo for std::logic_error ()
    from /lib64/libstdc++.so.6

(gdb) bt
#0  0x00007ffff56752c8 in typeinfo for std::logic_error ()
    from /lib64/libstdc++.so.6
#1  0x00007fffedd17492 in __cxa_decrement_exception_refcount ()
    from /usr/local/lib64/libc++abi.so.1
#2  0x00007fffedd6b24f in 
std::exception_ptr::operator=(std::exception_ptr const&) () from 
/usr/local/lib64/libc++.so.1
#3  0x00007fffe8ccd80d in SPLITT::ThreadExceptionHandler::CaptureException (
     this=0x1a2a78b8) at ./SPLITT.h:1853
#4 
SPLITT::ThreadExceptionHandler::Run<SPLITT::PostOrderTraversal<PCMBaseCpp::BM>::TraverseTreeSingleThreadLoopPrunes()::{lambda()#1}>(SPLITT::PostOrderTraversal<PCMBaseCpp::BM>::TraverseTreeSingleThreadLoopPrunes()::{lambda()#1}) 
(
     this=0x1a2a78b8, f=...) at ./SPLITT.h:1847
#5  0x00007fffe8cccf5d in 
SPLITT::PostOrderTraversal<PCMBaseCpp::BM>::TraverseTreeSingleThreadLoopPrunes 
(this=0x1a2a77f0) at ./SPLITT.h:2031
#6  0x00007fffe8cccafd in 
SPLITT::PostOrderTraversal<PCMBaseCpp::BM>::TraverseTree 
(this=0x1a2a77f0, mode=<optimized out>) at ./SPLITT.h:1884
#7  0x00007fffe8ccc9ef in 
SPLITT::TraversalTask<PCMBaseCpp::BM>::TraverseTree (
     this=0x1a2a5930, par=..., mode=4294662200) at ./SPLITT.h:582
#8  0x00007fffe8c47d67 in 
PCMBaseCpp::TraversalTaskWrapper<PCMBaseCpp::BM>::TraverseTree 
(this=0x7ffffffb5838, par=..., mode=0) at ./QuadraticPolyCommon.h:127
#9  0x00007fffe8cce193 in 
Rcpp::CppMethod2<PCMBaseCpp::TraversalTaskWrapper<PCMBaseCpp::BM>, 
std::__1::basic_string<char, std::__1::char_traits<char>, 
std::__1::allocator<char> >, std::__1::vector<double, 
std::__1::allocator<double> > const&, unsigned int>::operator() 
(this=0x9197280, object=<optimized out>,
     args=<optimized out>)
     at 
/data/gannet/ripley/R/test-clang/Rcpp/include/Rcpp/module/Module_generated_CppMethod.h:195
#10 0x00007fffe8cc9309 in 
Rcpp::class_<PCMBaseCpp::TraversalTaskWrapper<PCMBaseCpp::BM> 
 >::invoke_notvoid (this=<optimized out>, method_xp=<optimized out>,
     object=0x1079b8b8, args=0x7ffffffb5ab0, nargs=2)
     at 
/data/gannet/ripley/R/test-clang/Rcpp/include/Rcpp/module/class.h:234
#11 0x00007fffe9bfcab5 in CppMethod__invoke_notvoid (args=<optimized out>)
     at module.cpp:220
#12 0x000000000049dd6b in do_External (call=0xb987510, op=0x889268,
     args=0xf33cac8, env=0xf33cb00)
     at /data/gannet/ripley/R/svn/R-devel/src/main/dotcode.c:572

Please correct ASAP and before Nov 28 to safely retain the package on CRAN.

-- 
Brian D. Ripley,                  ripley@stats.ox.ac.uk
Emeritus Professor of Applied Statistics, University of Oxford


-- 
Venelin MITOV


-- 
Venelin MITOV