/* Form definition file generated with fdesign. */

#include "forms.h"
#include <stdlib.h>
#include "ggui.h"

FD_ggui *create_form_ggui(void)
{
  FL_OBJECT *obj;
  FD_ggui *fdui = (FD_ggui *) fl_calloc(1, sizeof(*fdui));

  fdui->ggui = fl_bgn_form(FL_NO_BOX, 770, 640);
  obj = fl_add_box(FL_UP_BOX,0,0,770,640,"");
  fdui->Browser = obj = fl_add_browser(FL_NORMAL_BROWSER,30,90,710,530,"");
    fl_set_object_lsize(obj,FL_NORMAL_SIZE);
    fl_set_object_lstyle(obj,FL_FIXED_STYLE);
  fdui->QuitButton = obj = fl_add_button(FL_NORMAL_BUTTON,600,30,70,40,"Quit");
    fl_set_object_lsize(obj,FL_MEDIUM_SIZE);
    fl_set_object_lstyle(obj,FL_BOLD_STYLE);
    fl_set_object_callback(obj,QuitButtonCB,0);
  fdui->StepButton = obj = fl_add_button(FL_NORMAL_BUTTON,90,30,70,40,"Step");
    fl_set_object_lsize(obj,FL_MEDIUM_SIZE);
    fl_set_object_lstyle(obj,FL_BOLD_STYLE);
    fl_set_object_callback(obj,StepButtonCB,0);
  fl_end_form();

  return fdui;
}
/*---------------------------------------*/

