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
