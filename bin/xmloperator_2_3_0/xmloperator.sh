#!/bin/bash
ulimit -s 2048
pushd ${XMLOPERATOR_USER_HOME}
java -cp ${XMLOPERATOR_HOME}/xmloperator.jar:${XALAN_LIB}/xml-apis.jar:${XALAN_LIB}/xercesImpl.jar:${XALAN_LIB}/xalan.jar org.xmloperator.Tool ${XMLOPERATOR_USER_HOME}/data/xmloperator.xin
popd
