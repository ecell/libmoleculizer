<?xml version="1.0"?>

<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
<xsl:output method="html"
            encoding="ISO-8859-1"/>

<!-- Generally convert species-tag into species-name. -->
<xsl:template match="species-tag">
 <xsl:apply-templates/>
</xsl:template>

<xsl:template match="species-tag" mode="itemMode">
<li><xsl:apply-templates/></li>
</xsl:template>

<!-- Convert plex-species into species-entry. -->
<xsl:template match="plex-species">
  <li><xsl:apply-templates select="./species-tag"/></li>
</xsl:template>

<!-- Convert stoch-species into species-entry. -->
<xsl:template match="stoch-species">
 <li><xsl:apply-templates select="./species-tag"/></li>
</xsl:template>

<xsl:template match="multiplicity">
multiplicity <xsl:apply-templates/>
</xsl:template>

<xsl:template match="reaction-substrate">
<li>substrate <xsl:apply-templates/> </li>
</xsl:template>

<xsl:template match="reaction-delta">
<li>delta <xsl:apply-templates/> </li>
</xsl:template>

<xsl:template match="reaction-rate">
<li>rate <xsl:apply-templates/> </li>
</xsl:template>

<xsl:template match="reaction">
<li> reaction <ul> <xsl:apply-templates/></ul></li>
</xsl:template>

<xsl:template match="dumpable">
 <li>dumpable <xsl:apply-templates select="./dumpable-name"/>
  <ul>
   <xsl:apply-templates select="./species-tag" mode="itemMode"/>
  </ul>
 </li>
</xsl:template>

<xsl:template match="mzr-inventory">
<HTML>
<HEAD><TITLE>Inventory Summary</TITLE></HEAD>
<body>
<p><H1>Inventory Summary</H1></p>
<p> <H2>Species</H2>
<ul>
 <xsl:apply-templates
  select="./unit-inventory/plex-inventory/plex-species"/>
 <xsl:apply-templates
  select="./unit-inventory/stoch-inventory/stoch-species"/>
</ul>
</p>
<p> <H2>Reactions</H2>
<ul>
 <xsl:apply-templates select="./reaction-inventory"/>
</ul>
</p>

<p> <H2>Dumpables</H2>
<ul>
 <xsl:apply-templates select="./dumpable-inventory"/>
</ul>
</p>
</body>
</HTML>
</xsl:template>

</xsl:stylesheet>



