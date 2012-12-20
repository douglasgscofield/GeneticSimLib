/* Form definition file generated with fdesign. */

#include "forms.h"
#include <stdlib.h>
#include "rgui.h"

FD_rgui *create_form_rgui(void)
{
  FL_OBJECT *obj;
  FD_rgui *fdui = (FD_rgui *) fl_calloc(1, sizeof(*fdui));

  fdui->rgui = fl_bgn_form(FL_NO_BOX, 670, 510);
  obj = fl_add_box(FL_UP_BOX,0,0,670,510,"");
  fdui->xyplot = obj = fl_add_xyplot(FL_FILL_XYPLOT,40,160,590,320,"");
    fl_set_object_boxtype(obj,FL_DOWN_BOX);
    fl_set_object_color(obj,FL_PALEGREEN,FL_RIGHT_BCOL);
    fl_set_object_lsize(obj,FL_MEDIUM_SIZE);
  fdui->AEntry = obj = fl_add_input(FL_FLOAT_INPUT,280,30,80,30,"A");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
    fl_set_object_lstyle(obj,FL_FIXED_STYLE);
  fdui->BEntry = obj = fl_add_input(FL_FLOAT_INPUT,280,70,80,30,"B");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
    fl_set_object_lstyle(obj,FL_FIXED_STYLE);
  fdui->NEntry = obj = fl_add_input(FL_INT_INPUT,280,110,80,30,"N");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
    fl_set_object_lstyle(obj,FL_FIXED_STYLE);
  fdui->DistributionChoice = obj = fl_add_choice(FL_NORMAL_CHOICE,100,30,100,30,"Distribution");
    fl_set_object_boxtype(obj,FL_SHADOW_BOX);
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
    fl_set_object_lstyle(obj,FL_BOLD_STYLE);
    fl_set_object_callback(obj,DistributionChoiceCB,0);
  fdui->PlotButton = obj = fl_add_button(FL_NORMAL_BUTTON,120,80,60,50,"Plot");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
    fl_set_object_lstyle(obj,FL_BOLD_STYLE);
    fl_set_object_callback(obj,PlotButtonCB,0);
  fdui->QuitButton = obj = fl_add_button(FL_NORMAL_BUTTON,580,30,50,50,"Quit");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
    fl_set_object_lstyle(obj,FL_BOLD_STYLE);
    fl_set_object_callback(obj,QuitButtonCB,0);
  fdui->ProgressSlider = obj = fl_add_slider(FL_HOR_FILL_SLIDER,60,180,550,20,"");
    fl_set_object_color(obj,FL_COL1,FL_BLUE);
     fl_set_slider_return(obj, FL_RETURN_CHANGED);
  fdui->MinEntry = obj = fl_add_input(FL_FLOAT_INPUT,420,70,80,30,"Min");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
    fl_set_object_lstyle(obj,FL_FIXED_STYLE);
  fdui->MaxEntry = obj = fl_add_input(FL_FLOAT_INPUT,420,110,80,30,"Max");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
    fl_set_object_lstyle(obj,FL_FIXED_STYLE);
  obj = fl_add_text(FL_NORMAL_TEXT,420,40,80,30,"Plot Range");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
    fl_set_object_lalign(obj,FL_ALIGN_CENTER|FL_ALIGN_INSIDE);
    fl_set_object_lstyle(obj,FL_BOLD_STYLE);
  fl_end_form();

  return fdui;
}
/*---------------------------------------*/

