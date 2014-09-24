smrtmake
========

Hackable smrtpipe workflows using makefiles instead of smrtpipe.py.  Mostly 
self contained, they offer benefits over smrpipe in several respects:

1. Restartable.  If something fails in the middle, you can (in many cases) restart where you left off. 
2. No XML.  The configuration and workflow is all contained in the Makefile.
3. Piecemeal execution.  Only execute the parts of the workflow you care about.
4. More freedom to organize outputs.
5. Customizable to your needs. Take out stuff you don't want or add stuff that you do want. 
Just write adapters to manage the interfacing.

Some downsides.

1. Make syntax is not exactly user friendly and you need to understand how command line works.
2. Command line only, no integration with SMRT Portal.
3. Limited support.

**NOTE**: These have only been tested on SGE clusters. You can probably find ways to run 
on other cluster setups with a little modification.

Quick Links
-----------

[Resequencing] (#reseq)  
[HGAP3] (#hgap3)
[Circularize] (#circ)

Versions
--------

As many of you know, we evolve the secondary software from time to time.  The makefiles 
will be branched accordingly so you can find makefiles compatible with the version of 
secondary that you have installed.  The master branch will be compatible with the most 
recently released version of secondary, usually taking advantage of some recent 
improvements to the code base that won't be found in older releases.

Running
-------

You'll need a smrtanalysis installation and SGE cluster manager in order to run 
the Makefile workflows, although it would be trivial to modify them to run on a 
single compute node.

All you need is two files: Makefile and input.fofn.

    > ls
    Makefile input.fofn
    > make

That's it!

It can also be parallelized using the '-j' option to *make*

    > make -j 8

Some people may just want to keep the makefiles in a central location (i.e., treat 
it like code).  In that case you'll have to supply the '-f' option to make and provide 
the full path to the make file.

    > make -f </path/to/file.mk>

If the makefile includes another makefile (e.g., reseq.mk includes chunk.mk), then 
you'll also need to include the '-I' option to specify the directory where it resides. 
Using the previous example:

    > make -f </path/to/file.mk> -I </directory/containing/included/makefiles>

The reseq.mk file falls into this category.  So if you just cloned the repository and 
wanted to 

**NOTE:** There are some assumptions builtin to this workflow that may not hold 
in your environment, e.g., metadata.xml file locations.  While fixable, it's not 
always obvious what you need to do.  Send me an email and I can try to help you 
get unstuck.

### chunk.mk

A simple recipe for splitting the input.fofn file into chunks that can be maintained 
in other makefiles.  This is included in other makefiles here, so you'll probably want 
to include this.

###<a name="reseq"/> reseq.mk

Runs the resequencing workflow, which consists of filter, mapping and quiver consensus. 
You'll need to specify the location of a reference and a File-Of-Filenames (FOFN) at 
minimum. This make file includes chunk.mk, so you'll need to have that in your path (see 
the section on 'Running')

###<a name="circ"/> circularize.mk

Given a polished assembly file (uncompressed), identify those that appear 
circular. 

    > ls
      circularize.mk 	polished_assembly.fasta
    > make -sf circularize.mk
      unitig_2: looks circular
      unitig_4: looks circular
      unitig_1: looks circular

You can also override the following defaults from the command line:
	FASTA: Location of the polished assembly (default ./polished_assembly.fasta)
	OVLSIZE: Determines how much sequence is used to detect an overlap.
 	         Sequences < 2 x OVLSIZE are ignored (default 20000)

Example:
    > make -sf circularize.mk FASTA=/path/to/polished_assembly.fasta

NOTE: This only detects perfectly circular assemblies.  Some imperfect 
assemblies that may be circular will likely be missed.

### <a name="hgap3"/> hgap3.mk

HGAP.3 workflow tuned for larger genomes.  It doesn't generate reports (though they may 
be added in the future) and the output is organized much differently, but better IMHO, 
than smrtpipe.

This was built and run using PacBio's internal cluster, which uses SGE.  The nodes were 
fairly sizeable with 48 CPU, 200GB RAM. Consider modifying the Makefile to suite your needs:

    > head hgap3.mk
    # SGE queue name to submit jobs to
    QUEUE ?= huasm
    # Estimated size of the genome to assemble
    GENOME_SIZE ?= 700000000
    # Splits data into this many chunks, each chunk processed independently
    CHUNK_SIZE ?= 15
    # How many threads a process will use (also how many SGE slots will be requested)
    NPROC ?= 32
    # Local temp root directory, must have write access and a decent amount of space (~100GB)
    LOCALTMP ?= /scratch

The default recipe will run the HGAP.3 workflow to completion, i.e.:

    > make -j 15

However, you may only want to run it to a certain step.  You can also continue where you left off.

    # runs the workflow to get corrected reads
    > make -j 15 correction

    # runs (or continues) the workflow to get a draft assembly
    > make draft

Dry-run is also supported, when you want to double check what will be run.

    # displays the list of commands that will be run in order
    > make -n

