#ifndef FD_rgui_h_
#define FD_rgui_h_
/* Header file generated with fdesign. */

/**** Callback routines ****/

extern void DistributionChoiceCB(FL_OBJECT *, long);
extern void PlotButtonCB(FL_OBJECT *, long);
extern void QuitButtonCB(FL_OBJECT *, long);


/**** Forms and Objects ****/

typedef struct {
	FL_FORM *rgui;
	FL_OBJECT *xyplot;
	FL_OBJECT *AEntry;
	FL_OBJECT *BEntry;
	FL_OBJECT *NEntry;
	FL_OBJECT *DistributionChoice;
	FL_OBJECT *PlotButton;
	FL_OBJECT *QuitButton;
	FL_OBJECT *ProgressSlider;
	FL_OBJECT *MinEntry;
	FL_OBJECT *MaxEntry;
	void *vdata;
	long ldata;
} FD_rgui;

extern FD_rgui * create_form_rgui(void);

#endif /* FD_rgui_h_ */
