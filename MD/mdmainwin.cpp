#include "mdmainwin.h"
#include "ui_mdmainwin.h"

mdmainwin::mdmainwin(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::mdmainwin)
{
    ui->setupUi(this);
}

mdmainwin::~mdmainwin()
{
    delete ui;
}

void mdmainwin::on_actionExit_triggered()
{
    QCoreApplication::exit();
}
