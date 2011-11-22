#include <QtGui/QApplication>
#include "mdmainwin.h"

int main(int argc, char *argv[])
{
    QApplication application(argc, argv);
    mdmainwin window;
    window.show();

    int result = application.exec();
    return result;
}
