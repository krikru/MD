#include <QApplication>
#include <QPushButton>

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    QPushButton hello("Hello World!");
    hello.resize(100, 30);
    QObject::connect(&hello, SIGNAL(clicked()),
        &app, SLOT(quit()));
    hello.show();
    return app.exec();
}