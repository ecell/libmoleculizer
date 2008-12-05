package integration;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import org.xmloperator.util.exception.ExceptionControler;
import org.xmloperator.framework.control.ExternalControl;
import org.xmloperator.Tool;

// This xmloperator integration example opens two documents
// within two instances of the tool.
// A "Close all" menu item added to the "Options" menu allows some work
// on these documents before closing the tool instances.
public class MyProgram2 {
  private org.xmloperator.framework.control.ExternalControl tool1, tool2;
  private final ActionListener closeAllListener =
    new ActionListener() {
      public void actionPerformed(ActionEvent event) {
        closeAll_actionPerformed();
      }
    };
  private final ActionListener instanceDeletedNotificationListener =
    new ActionListener() {
      public void actionPerformed(ActionEvent event) {
        notifyInstanceDeleted();
      }
    };

  public static void main(String[] args) {
    org.xmloperator.util.exception.ExceptionControler.software = "MyProgram";
    try {
      new MyProgram2();
    }
    catch(Exception exception) {
      org.xmloperator.util.exception.ExceptionControler.printStackTrace(
        exception);
      System.exit(1);
    }
  }

  private MyProgram2() {
    String myDocumentFile1 = "data/xmloperator_doc.xml";
    String myDocumentFile2 = "data/xmloperator_doc_fr.xml";

    tool1 = org.xmloperator.Tool.createToolInstance(myDocumentFile1);
    tool1.addMenuItem("Close All", closeAllListener, null);
    tool1.setInstanceDeletedNotificationListener(
      instanceDeletedNotificationListener);

    tool2 = org.xmloperator.Tool.createToolInstance(myDocumentFile2);
    tool2.addMenuItem("Close All", closeAllListener, null);
    tool2.setInstanceDeletedNotificationListener(
      instanceDeletedNotificationListener);
  }

  private void closeAll_actionPerformed() {
    try {
      //
      // Some work on tool1 and tool2 documents ...
      //
      tool1.deleteInstance();
      tool2.deleteInstance();
    }
    catch(Exception exception) {
      org.xmloperator.util.exception.ExceptionControler.printStackTrace(
        exception);
    }
  }

  private void notifyInstanceDeleted() {
    if (tool1.getInstanceCount() == 0)
      System.exit(0);
  }
}