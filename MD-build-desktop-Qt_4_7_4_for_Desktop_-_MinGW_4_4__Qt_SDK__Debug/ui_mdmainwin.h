/********************************************************************************
** Form generated from reading UI file 'mdmainwin.ui'
**
** Created: Thu 10. Nov 13:56:15 2011
**      by: Qt User Interface Compiler version 4.7.4
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MDMAINWIN_H
#define UI_MDMAINWIN_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QHeaderView>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QPushButton>
#include <QtGui/QRadioButton>
#include <QtGui/QStatusBar>
#include <QtGui/QToolBar>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_mdmainwin
{
public:
    QAction *actionExit;
    QAction *actionAbout;
    QWidget *centralWidget;
    QPushButton *pushButton;
    QRadioButton *radioButton;
    QMenuBar *menuBar;
    QMenu *menuFile;
    QMenu *menuHelp;
    QMenu *menuSettings;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *mdmainwin)
    {
        if (mdmainwin->objectName().isEmpty())
            mdmainwin->setObjectName(QString::fromUtf8("mdmainwin"));
        mdmainwin->resize(388, 278);
        actionExit = new QAction(mdmainwin);
        actionExit->setObjectName(QString::fromUtf8("actionExit"));
        actionAbout = new QAction(mdmainwin);
        actionAbout->setObjectName(QString::fromUtf8("actionAbout"));
        centralWidget = new QWidget(mdmainwin);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        pushButton = new QPushButton(centralWidget);
        pushButton->setObjectName(QString::fromUtf8("pushButton"));
        pushButton->setGeometry(QRect(30, 120, 75, 23));
        radioButton = new QRadioButton(centralWidget);
        radioButton->setObjectName(QString::fromUtf8("radioButton"));
        radioButton->setGeometry(QRect(240, 120, 82, 18));
        mdmainwin->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(mdmainwin);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 388, 18));
        menuFile = new QMenu(menuBar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        menuHelp = new QMenu(menuBar);
        menuHelp->setObjectName(QString::fromUtf8("menuHelp"));
        menuSettings = new QMenu(menuBar);
        menuSettings->setObjectName(QString::fromUtf8("menuSettings"));
        mdmainwin->setMenuBar(menuBar);
        mainToolBar = new QToolBar(mdmainwin);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        mdmainwin->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(mdmainwin);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        mdmainwin->setStatusBar(statusBar);

        menuBar->addAction(menuFile->menuAction());
        menuBar->addAction(menuSettings->menuAction());
        menuBar->addAction(menuHelp->menuAction());
        menuFile->addSeparator();
        menuFile->addAction(actionExit);
        menuHelp->addAction(actionAbout);

        retranslateUi(mdmainwin);

        QMetaObject::connectSlotsByName(mdmainwin);
    } // setupUi

    void retranslateUi(QMainWindow *mdmainwin)
    {
        mdmainwin->setWindowTitle(QApplication::translate("mdmainwin", "mdmainwin", 0, QApplication::UnicodeUTF8));
        actionExit->setText(QApplication::translate("mdmainwin", "Exit", 0, QApplication::UnicodeUTF8));
        actionAbout->setText(QApplication::translate("mdmainwin", "About", 0, QApplication::UnicodeUTF8));
        pushButton->setText(QApplication::translate("mdmainwin", "PushButton", 0, QApplication::UnicodeUTF8));
        radioButton->setText(QApplication::translate("mdmainwin", "RadioButton", 0, QApplication::UnicodeUTF8));
        menuFile->setTitle(QApplication::translate("mdmainwin", "File", 0, QApplication::UnicodeUTF8));
        menuHelp->setTitle(QApplication::translate("mdmainwin", "Help", 0, QApplication::UnicodeUTF8));
        menuSettings->setTitle(QApplication::translate("mdmainwin", "Settings", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class mdmainwin: public Ui_mdmainwin {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MDMAINWIN_H
