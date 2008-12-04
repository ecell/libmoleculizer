package integration;

import org.xmloperator.Tool;

// An example to integrate xmloperator within a Java program
public class MyProgram1 {

  public static void main(String[] args) {
    String myDocumentFile = "data/xmloperator_doc.xml";
    org.xmloperator.Tool.main(new String[] {myDocumentFile});
  }
}