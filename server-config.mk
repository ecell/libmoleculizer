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

# This is used only in server-based delivery.  It is used to set up the
# user's connection to "live" documentation in xmloperator, and it is used
# in server scripts.
#
# Users that will only be working with the standalone version should
# not have to worry about setting this.  They will get live documentation
# from static pages in the delivery directory.

MOLECULIZER_SERVER := http://genome.molsci.org

# Apache's main directory for web pages.  This is used to install a link to
# Moleculizer web pages and is used in server scripts.
#
# Note that Apache must be configured to follow links.
#
# Users that will only be working with the standalone version should not
# have to worry about setting this.

SERVER_HTDOCS_DIR := /srv/www/htdocs

# Apache's main directory for cgi scripts.
# This is used to install a cgi dispatch script and is used in server scripts.
#
# Users that will only be working with the standalone version should not
# have to worry about setting this.

SERVER_CGI_DIR := /srv/www/cgi-bin



