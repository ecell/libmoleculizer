<?xml version="1.0"?><!--*- XML -*-->

<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
<!-- $Id: RngToRncText.xsl,v 1.1 2003/02/07 12:30:25 daro Exp $ -->
<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->

<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     Copyright (c) 2002, Pantor Engineering AB
     All rights reserved.
     
     Redistribution and use in source and binary forms, with or
     without modification, are permitted provided that the following
     conditions are met:
     
     * Redistributions of source code must retain the above copyright
       notice, this list of conditions and the following disclaimer. 

     * Redistributions in binary form must reproduce the above
       copyright notice, this list of conditions and the following
       disclaimer in the documentation and/or other materials provided
       with the distribution.

     * Neither the name of Pantor Engineering AB nor the names of its
       contributors may be used to endorse or promote products derived
       from this software without specific prior written permission.
     
     THIS   SOFTWARE   IS PROVIDED BY    THE   COPYRIGHT  HOLDERS  AND
     CONTRIBUTORS "AS IS"   AND ANY  EXPRESS OR  IMPLIED   WARRANTIES,
     INCLUDING,  BUT  NOT LIMITED  TO,   THE  IMPLIED  WARRANTIES   OF
     MERCHANTABILITY    AND  FITNESS  FOR   A  PARTICULAR  PURPOSE ARE
     DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
     BE  LIABLE   FOR ANY    DIRECT, INDIRECT,   INCIDENTAL,  SPECIAL,
     EXEMPLARY, OR CONSEQUENTIAL DAMAGES  (INCLUDING, BUT NOT  LIMITED
     TO, PROCUREMENT  OF  SUBSTITUTE GOODS OR  SERVICES;  LOSS OF USE,
     DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
     ANY THEORY OF  LIABILITY, WHETHER IN CONTRACT,  STRICT LIABILITY,
     OR  TORT (INCLUDING NEGLIGENCE OR  OTHERWISE) ARISING  IN ANY WAY
     OUT OF  THE  USE OF   THIS  SOFTWARE,  EVEN IF ADVISED   OF   THE
     POSSIBILITY OF SUCH DAMAGE.

     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->

<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     Created by David.Rosenborg@pantor.com
     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->

<!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     RngToRncText.xsl converts a RELAX NG schema in XML syntax to the
     compact syntax.

     If the source document is a RELAX NG schema, an extension
     function converting result tree fragments to node sets is
     required. If the processor supports the exslt:node-set function,
     no configuration should be needed. See RngToRncProcessorConfig.xsl
     for further information.

     This stylesheet can also be applied to the intermediate XML
     representation of a compact schema as output by the
     RngToRncXml.xsl stylesheet. In this case, no extension function
     is required.

     For a description of the underlying XML to compact syntax
     transformation, see RngToRncXml.xsl.

     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->

<xsl:transform
  version="1.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:rng="http://relaxng.org/ns/structure/1.0"
  xmlns:pax="http://www.pantor.com/pax"
>

  <xsl:import href="RngToRncXml.xsl"/>
  <xsl:import href="RngToRncProcessorConfig.xsl"/>

  <xsl:output method="text"/>

  <!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->
  <!-- Parameters -->
  <!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->

  <!-- collapse-lines:
     If true, output constructs spanning multiple lines will be
     groupd into a single line unless it exceeds $line-width chars. -->

  <xsl:param name="collapse-lines" select="true ()"/>

  <!-- indent-width:
     The number of characters to indent at each indentation level -->

  <xsl:param name="indent-width" select="3"/>

  <!-- line-width:
     see the group-lines parameter. -->

  <xsl:param name="line-width" select="80"/>

  <!-- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ -->

  <xsl:template match="/">
    <xsl:if test="$recursive">
      <xsl:message terminate="yes">
	<xsl:text>ERROR: recursive processing is not </xsl:text>
	<xsl:text>possible with the text output method</xsl:text>
      </xsl:message>
    </xsl:if>

    <xsl:choose>
      <xsl:when test="rng:*">
	<xsl:call-template name="make-body-from-r-t-f">
	  <xsl:with-param name="schema">
	    <xsl:apply-imports/>
	  </xsl:with-param>
	</xsl:call-template>
      </xsl:when>

      <xsl:otherwise>
	<xsl:call-template name="make-body"/>
      </xsl:otherwise>
    </xsl:choose>

    <xsl:text>&#10;</xsl:text>
  </xsl:template>

  <xsl:template name="make-body">
    <xsl:apply-templates mode="keep"/>
  </xsl:template>

  <xsl:template match="group" mode="keep">
    <xsl:choose>
      <xsl:when 
	test="
	not ($collapse-lines) or 
	@collapse = 'no' or
	.//doc or 
	.//group [@collapse = 'no']">
	<xsl:apply-templates mode="keep"/>
      </xsl:when>
      <xsl:otherwise>
	<xsl:variable name="size" select="sum (.//*/@size)"/>
	<xsl:variable name="level" select="count (ancestor::indent)"/>
	<xsl:choose>
	  <xsl:when test="$line-width > ($size + ($level * $indent-width))">
	    <xsl:apply-templates mode="flatten"/>
	  </xsl:when>
	  <xsl:otherwise>
	    <xsl:apply-templates mode="keep"/>
	  </xsl:otherwise>
	</xsl:choose>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <xsl:template match="sp" mode="keep">
    <xsl:text> </xsl:text>
  </xsl:template>

  <xsl:variable name="spaces" 
    select="concat (
    '                                        ',
    '                                        '
    )"/>

  <xsl:template match="nl" mode="keep">
    <xsl:text>&#10;</xsl:text>
    <xsl:variable name="level" select="count (ancestor::indent)"/>
    <xsl:variable name="following-op"
      select="following-sibling::*[1][self::op]"/>
    <xsl:choose>
      <xsl:when test="$following-op">
	<xsl:value-of
	  select="substring ($spaces, 1, 
	  $level * $indent-width - $following-op/@size)"/>
      </xsl:when>
      <xsl:otherwise>
	<xsl:value-of 
	  select="substring ($spaces, 1, $level * $indent-width)"/>
      </xsl:otherwise>
    </xsl:choose>
	
  </xsl:template>

  <xsl:template match="sp | nl" mode="flatten">
    <xsl:text> </xsl:text>
  </xsl:template>
</xsl:transform>
