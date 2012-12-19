#ifndef FD_ggui_h_
#define FD_ggui_h_
/* Header file generated with fdesign. */

/**** Callback routines ****/

extern void QuitButtonCB(FL_OBJECT *, long);
extern void StepButtonCB(FL_OBJECT *, long);


/**** Forms and Objects ****/

typedef struct {
	FL_FORM *ggui;
	FL_OBJECT *Browser;
	FL_OBJECT *QuitButton;
	FL_OBJECT *StepButton;
	void *vdata;
	long ldata;
} FD_ggui;

extern FD_ggui * create_form_ggui(void);

#endif /* FD_ggui_h_ */
