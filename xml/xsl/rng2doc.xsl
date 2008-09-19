<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns="http://relaxng.org/ns/structure/1.0" xmlns:rng="http://relaxng.org/ns/structure/1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
><xsl:output method="xml"
/><xsl:template match="node() | @*"
><xsl:copy
><xsl:apply-templates select="node() | @*"
/></xsl:copy></xsl:template><xsl:template match="rng:grammar"
><xsl:element name="grammar-doc"
><xsl:apply-templates
/></xsl:element></xsl:template></xsl:stylesheet>