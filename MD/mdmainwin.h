#ifndef MDMAINWIN_H
#define MDMAINWIN_H

//Standard includes
//#include <vector>
//#include <string>

// Own includes
#include "mdsystem.h"
#include "settings.h"

// Qt includes
#include <QMainWindow>

namespace Ui {
    class mdmainwin;
}

class mdmainwin : public QMainWindow
{
    Q_OBJECT

private:
    //Private variables
    settings system_settings;

public:
    // Constructor and destructor
    explicit mdmainwin(QWidget *parent = 0);
    ~mdmainwin();

private slots:

    void closeEvent(QCloseEvent *event);

    // Line edits
    void on_sigma_le_editingFinished               ();
    void on_epsilon_le_editingFinished             ();
    void on_mass_le_editingFinished                ();
    void on_lattice_constant_le_editingFinished    ();
    void on_num_particles_le_editingFinished       ();
    void on_init_temperature_le_editingFinished    ();
    void on_desired_pressure_le_editingFinished    ();
    void on_desire_temperature_le_editingFinished  ();
    void on_time_step_le_editingFinished           ();
    void on_num_time_steps_le_editingFinished      ();
    void on_inner_cutoff_le_editingFinished        ();
    void on_outer_cutoff_le_editingFinished        ();
    // Spin boxes
    void on_measurement_interval_sb_editingFinished();
    // Radio buttons
    void on_npe_rb_clicked();
    void on_nve_rb_clicked();
    void on_nvt_rb_clicked();
    void on_npt_rb_clicked();
    // Combo boxes
    void on_lattice_type_cb_activated         (const QString &arg1);
    void on_epsilon_unit_cb_activated         (const QString &arg1);
    void on_desired_pressure_unit_cb_activated(const QString &arg1);
    // Check boxes
    void on_diffoceff_cb_clicked                (bool checked);
    void on_pressure_cb_clicked                 (bool checked);
    void on_cv_cb_clicked                       (bool checked);
    void on_msd_cb_clicked                      (bool checked);
    void on_energy_total_cb_clicked             (bool checked);
    void on_energy_kinetic_cb_clicked           (bool checked);
    void on_energy_potential_cb_clicked         (bool checked);
    void on_cohesive_energy_cb_clicked          (bool checked);
    void on_store_particle_possitions_cb_clicked(bool checked);
    //void on_draw_particles_cb_clicked           (bool checked);
    // Push buttons
    void on_save_element_pb_clicked    ();
    void on_load_element_pb_clicked    ();
    void on_start_simulation_pb_clicked();
    // Button boxes
    void on_settings_bb_accepted();
    void on_settings_bb_rejected();

private:
    // Private functions
    void write_to_text_browser(string output);

    // Static private functions
    static void static_write_to_text_browser(void* void_ptr_mainwin, string output);
    static void static_process_events       (void* void_ptr_mainwin               );

    // Private variables
    mdsystem simulation;

private:
    Ui::mdmainwin *ui;

};

#endif // MDMAINWIN_H
