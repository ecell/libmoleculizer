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

# Display a particular xml file with (iconic) edit links on all editable
# attributes.

use CGI;

# This seems to be the way that I invoke the XSLT transformation.  I then have
# to cat the (temporary) output file from this transformation to stdout and
# delete the output file.

sub writeVariationForm
{
    my ($session,
	$variation_source,
	$variation_target,
	$schema,
	$caption,
	$proxy_target,
	$output_file_path) = @_;

    $path_doc_url
	= $cgi_dispatch_url
	. "?task="
	. $mzr_doc_query_task
	. "&"
	. $schema_field
	. "="
	. $schema
	. "&"
	. $xpath_field
	. "=";

    $session_dir = "$sessions_dir/$session";
    $source_file_path = "$session_dir/$variation_source";

    # Run the XSLT transformation.
    system("java", "org.apache.xalan.xslt.Process",
	   "-in", "$source_file_path",
	   "-xsl", "$xsl_dir/mzr2form.xsl",
	   "-param", "action-url", "$cgi_dispatch_url",
	   "-param", "path-doc-url", "$path_doc_url",
	   "-param", "caption", "$caption",
	   "-param", "schema-doc-file-path", "$schema_doc_dir/$schema",
	   "-param", "variation-source", "$variation_source",
	   "-param", "variation-target", "$variation_target",
	   "-param", "session", "$session",
	   "-param", "proxy-target", "$proxy_target",
	   "-html",
	   "-out", "$output_file_path");
}


# This seems to be the place in the main body of the old cgi where I invoke
# the above routine to do the XSLT transformation.

elsif($task eq $variation_form_task)
{
    $session = $query->param($session_field);
    $session_dir = "$sessions_dir/$session";

    $variation_source = $query->param($variation_source_field);
    $schema = $query->param($schema_field);
    $variation_caption = $query->param($variation_caption_field);
    $variation_target = $query->param($variation_target_field);
    $proxy_target = $query->param($proxy_target_field);

    print STDERR "variation-form: proxy target is $proxy_target.\n";

    # Write variation form to temp file, because Xalan always puts its output
    # in a file.
    $form_tmp_file = "/tmp/variation-form-$$";
    writeVariationForm($session,
		       $variation_source,
		       $variation_target,
		       $schema,
		       $variation_caption,
		       $proxy_target,
		       $form_tmp_file);

    # Print the content-type header; defaulting to html.  Note that this
    # has to be done in each task that actually emits a web page (rather
    # than a redirection.)
    print $query->header();

    # Emit tmp file.
    catFile($form_tmp_file);

    # Remove tmp file.
    unlink($form_tmp_file);
}
