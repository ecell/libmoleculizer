<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns="http://www.w3.org/1999/xhtml" xmlns:htm="http://www.w3.org/1999/xhtml" xmlns:jsk="http://www.xmloperator.net/namespace/java" xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
><xsl:param name="title"
/><xsl:template match="/jsk:javaSkeleton"
><html xmlns="http://www.w3.org/1999/xhtml"
><head
><title
><xsl:value-of select="$title"
/></title></head><body
><h1
><xsl:value-of select="$title"
/></h1><div
><xsl:for-each select="jsk:import"
><xsl:call-template name="packageList"
><xsl:with-param name="namePrefix" select="&quot;import &quot;"
/><xsl:with-param name="idPrefix" select="&quot;_IMPORT_&quot;"
/></xsl:call-template></xsl:for-each><xsl:call-template name="packageList"
/><xsl:for-each select="jsk:import"
><xsl:call-template name="packages"
><xsl:with-param name="namePrefix" select="&quot;import &quot;"
/><xsl:with-param name="idPrefix" select="&quot;_IMPORT_&quot;"
/></xsl:call-template></xsl:for-each><xsl:call-template name="packages"
/></div><hr
/></body></html></xsl:template><xsl:template name="packageList"
><xsl:param name="namePrefix" select="&quot;&quot;"
/><xsl:param name="idPrefix" select="&quot;&quot;"
/><xsl:if test="jsk:package | jsk:defaultPackage"
><div
><ul
><xsl:for-each select="jsk:package"
><li
><xsl:value-of select="$namePrefix"
/><xsl:text
>package </xsl:text><a href="#{concat($idPrefix, @qualifiedIdentifier)}"
><i
><xsl:value-of select="@qualifiedIdentifier"
/></i></a></li></xsl:for-each><xsl:for-each select="jsk:defaultPackage"
><li
><a href="#{concat($idPrefix, &quot;_DEFAULT_&quot;)}"
><xsl:value-of select="$namePrefix"
/><xsl:text
>default package</xsl:text></a></li></xsl:for-each></ul></div></xsl:if></xsl:template><xsl:template name="packages"
><xsl:param name="namePrefix" select="&quot;&quot;"
/><xsl:param name="idPrefix" select="&quot;&quot;"
/><xsl:for-each select="jsk:package"
><xsl:call-template name="package"
><xsl:with-param name="namePrefix" select="$namePrefix"
/><xsl:with-param name="idPrefix" select="$idPrefix"
/><xsl:with-param name="id" select="@qualifiedIdentifier"
/></xsl:call-template></xsl:for-each><xsl:for-each select="jsk:defaultPackage"
><xsl:call-template name="package"
><xsl:with-param name="namePrefix" select="$namePrefix"
/><xsl:with-param name="idPrefix" select="$idPrefix"
/><xsl:with-param name="id" select="&quot;_DEFAULT_&quot;"
/></xsl:call-template></xsl:for-each></xsl:template><xsl:template name="package"
><xsl:param name="namePrefix"
/><xsl:param name="idPrefix"
/><xsl:param name="id"
/><xsl:variable name="pid" select="concat($idPrefix, $id)"
/><h3
><a id="{$pid}" name="{$pid}"
><xsl:value-of select="$namePrefix"
/><xsl:choose
><xsl:when test="$id = &quot;_DEFAULT_&quot;"
><xsl:text
>default package</xsl:text></xsl:when><xsl:otherwise
><xsl:text
>package </xsl:text><i
><xsl:value-of select="$id"
/></i></xsl:otherwise></xsl:choose></a></h3><xsl:call-template name="packageList"
><xsl:with-param name="namePrefix" select="$namePrefix"
/><xsl:with-param name="idPrefix" select="$idPrefix"
/></xsl:call-template><xsl:call-template name="classeList"
/><xsl:call-template name="uses"
/><xsl:call-template name="packages"
><xsl:with-param name="namePrefix" select="$namePrefix"
/><xsl:with-param name="idPrefix" select="$idPrefix"
/></xsl:call-template><xsl:call-template name="classes"
/></xsl:template><xsl:template name="classeList"
><xsl:if test="jsk:classOrInterface"
><div
><ul
><xsl:for-each select="jsk:classOrInterface"
><li
><xsl:text
>class or interface </xsl:text><a href="#{@qualifiedIdentifier}"
><i
><xsl:value-of select="@qualifiedIdentifier"
/></i></a></li></xsl:for-each></ul></div></xsl:if></xsl:template><xsl:template name="classes"
><xsl:for-each select="jsk:classOrInterface"
><xsl:call-template name="class"
/></xsl:for-each></xsl:template><xsl:template name="class"
><h4
><a id="{@qualifiedIdentifier}" name="{@qualifiedIdentifier}"
><xsl:text
>class or interface </xsl:text><i
><xsl:value-of select="@qualifiedIdentifier"
/></i></a></h4><xsl:call-template name="uses"
/></xsl:template><xsl:template name="uses"
><xsl:for-each select="jsk:uses"
><xsl:call-template name="use"
><xsl:with-param name="usage" select="&quot;uses&quot;"
/></xsl:call-template></xsl:for-each><xsl:for-each select="jsk:usedBy"
><xsl:call-template name="use"
><xsl:with-param name="usage" select="&quot;used by&quot;"
/></xsl:call-template></xsl:for-each></xsl:template><xsl:template name="use"
><xsl:param name="usage"
/><xsl:if test="jsk:type"
><div
><dl
><dd
><dl
><xsl:for-each select="jsk:type"
><dt
><xsl:value-of select="$usage"
/><xsl:text
> </xsl:text><a href="#{@qualifiedIdentifier}"
><xsl:value-of select="@qualifiedIdentifier"
/></a></dt></xsl:for-each></dl></dd></dl></div></xsl:if></xsl:template></xsl:stylesheet>