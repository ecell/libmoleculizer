#!/usr/bin/perl
###############################################################################
# Moleculizer - a stochastic simulator for cellular chemistry.
# Copyright (C) 2001  Walter Lawrence (Larry) Lok.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#    
# Contact information:
#   Larry Lok, Research Fellow          Voice: 510-981-8740
#   The Molecular Sciences Institute      Fax: 510-647-0699
#   2168 Shattuck Ave.                  Email: lok@molsci.org
#   Berkeley, CA 94704
###############################################################################

# This script makes targets for a simulation session's HTML menu.
# Which targets are made depends on the status of the Condor job
# associated to the session.

# For now, we're just interested in whether the job is running and whether
# it has completed.

# Subroutines to make the target files.
sub canViewJob
{
    system("touch can-view-job");
}

sub canDownloadJob
{
    system("touch can-download-job");
}

sub canListQueue
{
    system("touch can-list-queue");
}

$queue_dir = $ENV{'SERVER_HTDOCS_DIR'} . "/moleculizer/queue";

$input_session = $ARGV[0];

print STDERR "condor-session-status.pl: Examining Condor queue for $input_session.\n";

# Read and invert cluster to session-info mapping.
open(CLUSTER_MAP, "$queue_dir/cluster_map")
    or die("Could not open cluster-to-session file $queue_dir/cluster_map.\n");
while(<CLUSTER_MAP>)
{
    chop;

    # Note here I'm using % as the field separator in the cluster-to-session
    # file, because must admit whitespace in catch phrase.
    ($cluster, $session, $catch_phrase, $reply_email) = split /%/;

    if($session eq $input_session)
    {
	$input_cluster = $cluster;
    }
}
close CLUSTER_MAP;

# Did we find the cluster associated to the session?
if($input_cluster)
{
    # Open a pipe to read output from condor_q.
    open(CONDOR_Q, "condor_q $cluster|")
	or die("Could not open pipe from condor_q.\n");

    # Read down through header line.
    while(<CONDOR_Q>)
    {
	if(m/ID/)
	{
	    last;
	}
    }

    # Try to read the next line, which won't exist if the job is not
    # in the queue.  Something more complex than this will be needed for
    # jobs with multiple processes...
    #
    # Note here that the "angle bracket operator" does not automatically
    # stuff the next line into $_ except in a "while" condition!  Very weird.
    if($_ = <CONDOR_Q>)
    {
	chop;

	# Here "cluster" is what condor_q emits, namely cluster.process.
	($cluster, $owner, $submit_date, $submit_time,
	 $run_time, $status, $priority, $size, $command) = split;

	if($status =~ m/U/)
	{
	    print STDERR "Job unexpanded.\n";
	}
	elsif($status =~ m/H/)
	{
	    print STDERR "Job holding.\n";
	}
	elsif($status =~ m/R/)
	{
	    print STDERR "Job running.\n";
	    canListQueue();
	    canViewJob();
	}
	elsif($status =~ m/I/)
	{
	    print STDERR "Job idle.\n";
	    canListQueue();
	}
	elsif($status =~ m/C/)
	{
	    print STDERR "Job completed.\n";
	    canViewJob();
	    canDownloadJob();
	}
	elsif($status =~ m/X/)
	{
	    print STDERR "Job removed.\n";
	    canViewJob();
	    canDownloadJob();
	}
	else
	{
	    print STDERR "Job queued with status $status.\n";
	}
    }
    else
    {
	print STDERR "Job is history.\n";
	canViewJob();
	canDownloadJob();
    }

    close CONDOR_Q;
}
else
{
    print STDERR "Session is unmapped.\n";
}
