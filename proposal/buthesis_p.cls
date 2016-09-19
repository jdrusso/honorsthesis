%%% LaTeX class file to produce format for Bucknell University
%%% Honors Thesis.  Modification of standard LaTeX book class.
%%% This is very preliminary.
%%% Martin Ligare, Bucknell Physics Department
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{buthesis}[2001/04/28 Bucknell thesis class]
\LoadClass[12pt,onecolumn,oneside]{article}
\RequirePackage{ifthen}
\newboolean{@twoadv}
\setboolean{@twoadv}{false}
\DeclareOption{numbib}{\PassOptionsToPackage{numbers}{natbib}}
\DeclareOption{twoadv}{\setboolean{@twoadv}{true}}
\ProcessOptions\relax
%
\RequirePackage{bm}
\RequirePackage[square]{natbib}
% \RequirePackage[bf,hang,small]{caption2}
\RequirePackage{indentfirst}
\pagestyle{plain}
%\markright{}
\setlength{\parskip}{0.2in}
\setlength{\topmargin}{0.0in}
\setlength{\oddsidemargin}{0.5in}
\setlength{\evensidemargin}{0.5in}
\setlength{\textwidth}{6.0in}
\addtolength{\textheight}{0.4in}
%\setlength{\footskip}{1.0in}
\newcommand{\advisor}[1]{\newcommand{\advisorname}{#1}}
\newcommand{\advisorb}[1]{\newcommand{\advisornameb}{#1}}
\newcommand{\chair}[1]{\newcommand{\chairname}{#1}}
\newcommand{\department}[1]{\newcommand{\departmentname}{#1}}
\newcommand{\butitle}[1]{\newcommand{\titletext}{#1}}
\newcommand{\degree}[1]{\newcommand{\degreename}{#1}}
\newcommand{\acknowledgments}[1]{\thispagestyle{myheadings}\markright{}
                                 \setcounter{page}{2}
                                 \mbox{}
                                 \vspace{1.5in}

                                  \begin{center}
                                  {\large\bf Acknowledgments}
                                  \end{center}

                                  #1}
\renewcommand*{\bibname}{References}
%%% Remove following after spring 2005 - Not necessary because page number
%%% location fixed in chapter environment
\newcommand{\butableofcontents}{
                    \tableofcontents\thispagestyle{headings}\markright{}
                    \pagestyle{headings}}
\newcommand{\bulistoffigures}{\pagestyle{myheadings}\markright{}
                    \listoffigures\thispagestyle{myheadings}\markright{}
                    \pagestyle{headings}}
\newcommand{\bulistoftables}{\pagestyle{myheadings}\markright{}
                    \listoftables\thispagestyle{myheadings}\markright{}
                    \pagestyle{headings}}
\newcommand{\buchapter}[1]{\pagestyle{empty}
                    \chapter{#1}\thispagestyle{myheadings}\markright{}
                    \pagestyle{headings}}
\renewcommand{\maketitle}{\begin{titlepage}
\mbox{}
%\addtolength{\textheight}{1.0in}
\begin{center}
\renewcommand{\baselinestretch}{1.2}
\large
{\bf \MakeUppercase{\titletext}}

\renewcommand{\baselinestretch}{1.}
\normalsize
\vspace{0.1in}

by\\
\vspace{0.5in}
\@author\\
\vspace{0.5in}
A Proposal Submitted to the Honors Council\\
\vspace{0.1in}
For Honors in the Department of \departmentname \\
\vspace{0.1in}
\today
\end{center}


\ifthenelse{\boolean{@twoadv}}{
\vspace{0,5in}

Approved by: \hspace{0.2in}\underline{\hspace{2.5in}}\\
\mbox{\hspace{1.3in}}\advisorname\\
\mbox{\hspace{1.3in}}Thesis Advisor
\vspace{0.25in}

\mbox{\hspace{1.0in}}\underline{\hspace{2.5in}}\\
\mbox{\hspace{1.3in}}\advisornameb\\
\mbox{\hspace{1.3in}}Thesis Advisor
\vspace{0.25in}

\mbox{\hspace{1.0in}}\underline{\hspace{2.5in}}\\
\mbox{\hspace{1.3in}}\chairname\\
\mbox{\hspace{1.3in}}Chair, Department of \departmentname}
{\vspace{1.0in}

Approved: \hspace{0.2in}\underline{\hspace{2.5in}}\\
\mbox{\hspace{1.3in}}\advisorname \\
\mbox{\hspace{1.3in}}Thesis Advisor
\vspace{0.5in}

\mbox{\hspace{1.0in}}\underline{\hspace{2.5in}}\\
\mbox{\hspace{1.3in}}\chairname\\
\mbox{\hspace{1.3in}}Chair, Department of \departmentname}
\end{titlepage}}
\endinput
