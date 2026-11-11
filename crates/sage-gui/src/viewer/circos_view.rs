//! The "Crosslink plot (Circos)" window of the Annotated spectra tab: the
//! loaded results' crosslinks as a circular plot, every choice of
//! [`CircosStyle`] at hand, and a preview rendered from the very SVG that is
//! exported. Styles are saved under a name as templates, in the library of
//! this computer, for other results to look the same.

use super::circos::{
    self, Background, CircosData, CircosStyle, ColorBy, Highlight, LabelMode, LegendCorner, LinkKind, Palette, Pick, Preset, Template,
};
use super::figure;
use super::model::Results;
use crate::library::Library;
use eframe::egui::{self, Color32, RichText, TextEdit, Ui};
use std::path::Path;
use std::sync::Arc;

pub struct CircosView {
    pub open: bool,
    style: CircosStyle,
    /// The plot's data, the results it came from and the selection it was built for
    data: Option<(Arc<CircosData>, usize, String)>,
    /// The preview, and what it shows: style, data and width in pixels
    preview: Option<(egui::TextureHandle, CircosStyle, usize, u32)>,
    format: figure::Format,
    dpi: u32,
    preset: Preset,
    message: Option<(bool, String)>,
    /// The saved template picked in the Templates menu, by name
    template: String,
    /// Name typed for a new template
    template_name: String,
    /// Save and Update keep the proteins too
    template_proteins: bool,
    /// Update, Rename (the name typed) or Delete waiting for a click
    confirm: Option<Confirm>,
    template_message: Option<(bool, String)>,
    /// The default template has been applied (when the window first opened)
    started: bool,
    /// The fonts installed (read when the font menu first opens)
    fonts: Option<Vec<String>>,
    /// Text typed to narrow the font menu
    font_filter: String,
    /// The highlight being set up in the Highlights section
    adding: Adding,
}

/// What a new highlight picks
#[derive(Clone, Copy, Default, PartialEq)]
enum AddWhat {
    #[default]
    Link,
    Protein,
    Between,
    Kind,
    Region,
}

impl AddWhat {
    const ALL: [AddWhat; 5] = [AddWhat::Link, AddWhat::Protein, AddWhat::Between, AddWhat::Kind, AddWhat::Region];

    fn label(self) -> &'static str {
        match self {
            AddWhat::Link => "one crosslink",
            AddWhat::Protein => "a protein's crosslinks",
            AddWhat::Between => "between two proteins",
            AddWhat::Kind => "a kind of crosslink",
            AddWhat::Region => "a stretch of a protein",
        }
    }
}

/// A highlight being set up: what it picks, and the choices made for it
struct Adding {
    what: AddWhat,
    link: Option<((String, usize), (String, usize))>,
    /// Text typed to narrow the crosslink menu
    filter: String,
    a: String,
    b: String,
    kind: LinkKind,
    start: usize,
    end: usize,
}

impl Default for Adding {
    fn default() -> Self {
        Adding {
            what: AddWhat::Link,
            link: None,
            filter: String::new(),
            a: String::new(),
            b: String::new(),
            kind: LinkKind::Multivalent,
            start: 1,
            end: 100,
        }
    }
}

impl Adding {
    /// The highlight's pick, once everything it needs is chosen
    fn pick(&self) -> Option<Pick> {
        match self.what {
            AddWhat::Link => self.link.clone().map(|(a, b)| Pick::Link { a, b }),
            AddWhat::Protein => (!self.a.is_empty()).then(|| Pick::Protein { protein: self.a.clone() }),
            AddWhat::Between => (!self.a.is_empty() && !self.b.is_empty()).then(|| Pick::Between { a: self.a.clone(), b: self.b.clone() }),
            AddWhat::Kind => Some(Pick::Kind { kind: self.kind }),
            AddWhat::Region => (!self.a.is_empty() && self.start <= self.end)
                .then(|| Pick::Region { protein: self.a.clone(), start: self.start, end: self.end }),
        }
    }
}

/// A menu of the plot's proteins (`(accession, label)`), setting `value`
fn protein_menu(ui: &mut Ui, id: &str, proteins: &[(String, String)], value: &mut String) {
    let shown = proteins.iter().find(|(a, _)| a == value).map(|(_, l)| l.clone()).unwrap_or_else(|| "protein".into());
    egui::ComboBox::from_id_salt(id).width(110.0).selected_text(shown).show_ui(ui, |ui| {
        for (accession, label) in proteins {
            ui.selectable_value(value, accession.clone(), label.as_str()).on_hover_text(accession.as_str());
        }
    });
}

/// The Highlights section: crosslinks picked out in colours of their own
fn highlights_ui(ui: &mut Ui, s: &mut CircosStyle, add: &mut Adding, data: &CircosData) {
    ui.label(
        RichText::new(
            "Crosslinks picked out in a colour of their own, drawn wider and on top. A crosslink takes the \
             colour of the first highlight that picks it.",
        )
        .weak(),
    );
    let counts: Vec<usize> =
        (0..s.highlights.len()).map(|i| data.links.iter().filter(|l| circos::highlight_of(s, l) == Some(i)).count()).collect();
    let hints: Vec<String> = s.highlights.iter().map(|h| h.pick.describe(data, s)).collect();
    let mut change: Option<(usize, i32)> = None;
    let n = s.highlights.len();
    if n > 0 {
        egui::Grid::new("circos_highlights").num_columns(5).spacing([6.0, 3.0]).show(ui, |ui| {
            for (i, h) in s.highlights.iter_mut().enumerate() {
                ui.checkbox(&mut h.on, "").on_hover_text("Drawn (untick to keep it without drawing it)");
                color_edit(ui, &mut h.color);
                ui.add_sized([150.0, ui.spacing().interact_size.y], TextEdit::singleline(&mut h.label).hint_text(&hints[i]))
                    .on_hover_text(format!("{}: the legend's text (empty: this)", hints[i]));
                ui.label(RichText::new(counts[i].to_string()).weak()).on_hover_text("crosslinks it colours");
                ui.horizontal(|ui| {
                    if ui.add_enabled(i > 0, egui::Button::new("⏶").small()).on_hover_text("before the one above").clicked() {
                        change = Some((i, -1));
                    }
                    if ui.add_enabled(i + 1 < n, egui::Button::new("⏷").small()).on_hover_text("after the one below").clicked() {
                        change = Some((i, 1));
                    }
                    if ui.add(egui::Button::new("remove").small()).clicked() {
                        change = Some((i, 0));
                    }
                });
                ui.end_row();
            }
        });
    }
    match change {
        Some((i, 0)) => {
            s.highlights.remove(i);
        }
        Some((i, -1)) => s.highlights.swap(i, i - 1),
        Some((i, _)) => s.highlights.swap(i, i + 1),
        None => {}
    }

    // a new one
    let proteins: Vec<(String, String)> = circos::shown_proteins(data, s).iter().map(|(p, label, _)| (p.accession.clone(), label.clone())).collect();
    ui.horizontal(|ui| {
        ui.label("add");
        egui::ComboBox::from_id_salt("circos_add_what").width(170.0).selected_text(add.what.label()).show_ui(ui, |ui| {
            for w in AddWhat::ALL {
                ui.selectable_value(&mut add.what, w, w.label());
            }
        });
    });
    ui.horizontal_wrapped(|ui| {
        match add.what {
            AddWhat::Link => {
                // every crosslink, most CSMs first, found by typing
                let mut links: Vec<&circos::Link> = data.links.iter().collect();
                links.sort_by(|a, b| b.csms.cmp(&a.csms).then(a.a.cmp(&b.a)));
                let text = |l: &circos::Link| {
                    format!("{} · {} CSM{}", Pick::Link { a: l.a.clone(), b: l.b.clone() }.describe(data, s), l.csms, if l.csms == 1 { "" } else { "s" })
                };
                let chosen = add
                    .link
                    .as_ref()
                    .and_then(|(a, b)| links.iter().find(|l| &l.a == a && &l.b == b))
                    .map(|l| text(l))
                    .unwrap_or_else(|| "pick a crosslink".into());
                egui::ComboBox::from_id_salt("circos_add_link")
                    .width(290.0)
                    .height(360.0)
                    .close_behavior(egui::PopupCloseBehavior::CloseOnClickOutside)
                    .selected_text(chosen)
                    .show_ui(ui, |ui| {
                        ui.add(TextEdit::singleline(&mut add.filter).hint_text("search: protein, residue").desired_width(270.0));
                        let needle = add.filter.trim().to_lowercase();
                        let mut picked = None;
                        for l in &links {
                            let t = text(l);
                            if !needle.is_empty() && !t.to_lowercase().contains(&needle) {
                                continue;
                            }
                            if ui.selectable_label(add.link.as_ref() == Some(&(l.a.clone(), l.b.clone())), t).clicked() {
                                picked = Some((l.a.clone(), l.b.clone()));
                            }
                        }
                        if picked.is_some() {
                            add.link = picked;
                            ui.memory_mut(|m| m.close_popup());
                        }
                    });
            }
            AddWhat::Protein => protein_menu(ui, "circos_add_a", &proteins, &mut add.a),
            AddWhat::Between => {
                protein_menu(ui, "circos_add_a", &proteins, &mut add.a);
                ui.label("and");
                protein_menu(ui, "circos_add_b", &proteins, &mut add.b);
            }
            AddWhat::Kind => {
                egui::ComboBox::from_id_salt("circos_add_kind").selected_text(add.kind.label()).show_ui(ui, |ui| {
                    for k in LinkKind::ALL {
                        ui.selectable_value(&mut add.kind, k, k.label());
                    }
                });
            }
            AddWhat::Region => {
                protein_menu(ui, "circos_add_a", &proteins, &mut add.a);
                ui.label("residues");
                ui.add(egui::DragValue::new(&mut add.start).range(1..=100_000));
                ui.label("to");
                ui.add(egui::DragValue::new(&mut add.end).range(1..=100_000));
            }
        }
        let pick = add.pick();
        if ui.add_enabled(pick.is_some(), egui::Button::new("Add")).clicked() {
            if let Some(pick) = pick {
                let color = circos::HIGHLIGHT_COLORS[s.highlights.len() % circos::HIGHLIGHT_COLORS.len()];
                s.highlights.push(Highlight::new(pick, color));
            }
        }
    });
    ui.horizontal(|ui| {
        ui.label("width ×").on_hover_text("Highlighted crosslinks: their width times this");
        ui.add(egui::DragValue::new(&mut s.highlight_width).range(0.2..=10.0).speed(0.05).max_decimals(2));
    });
    ui.horizontal(|ui| {
        ui.checkbox(&mut s.dim_others, "dim the others").on_hover_text("The crosslinks no highlight picks, in this colour and fainter");
        color_edit(ui, &mut s.dim_color);
    });
}

#[derive(Clone)]
enum Confirm {
    Update,
    Rename(String),
    Delete,
}

fn outcome(result: Result<(), String>, done: String) -> Option<(bool, String)> {
    Some(match result {
        Ok(()) => (false, done),
        Err(e) => (true, e),
    })
}

fn show_message(ui: &mut Ui, message: &Option<(bool, String)>) {
    if let Some((error, text)) = message {
        if *error {
            ui.label(RichText::new(format!("✖ {text}")).color(Color32::from_rgb(210, 60, 60)));
        } else {
            ui.label(RichText::new(format!("✔ {text}")).weak());
        }
    }
}

/// A template's name in the menu, marked if it is the default
fn template_label(lib: &Library, i: usize) -> String {
    let name = &lib.circos_templates[i].name;
    if lib.is_default_template(i) {
        format!("★ {name}")
    } else {
        name.clone()
    }
}

impl Default for CircosView {
    fn default() -> Self {
        CircosView {
            open: false,
            style: CircosStyle::default(),
            data: None,
            preview: None,
            format: figure::Format::Svg,
            dpi: 300,
            preset: Preset::Default,
            message: None,
            template: String::new(),
            template_name: String::new(),
            template_proteins: false,
            confirm: None,
            template_message: None,
            started: false,
            fonts: None,
            font_filter: String::new(),
            adding: Adding::default(),
        }
    }
}

/// `#rrggbb` ↔ the colour picker
fn color_edit(ui: &mut Ui, hex: &mut String) -> bool {
    let (r, g, b) = circos::rgb(hex);
    let mut color = Color32::from_rgb(r, g, b);
    let changed = ui.color_edit_button_srgba(&mut color).changed();
    if changed {
        *hex = format!("#{:02x}{:02x}{:02x}", color.r(), color.g(), color.b());
    }
    changed
}

fn section(ui: &mut Ui, title: &str, open: bool, body: impl FnOnce(&mut Ui)) {
    egui::CollapsingHeader::new(RichText::new(title).strong()).default_open(open).show(ui, body);
}

fn number(ui: &mut Ui, label: &str, value: &mut f64, range: std::ops::RangeInclusive<f64>, speed: f64, tip: &str) {
    ui.horizontal(|ui| {
        ui.label(label).on_hover_text(tip);
        ui.add(egui::DragValue::new(value).range(range).speed(speed).max_decimals(3));
    });
}

impl CircosView {
    /// The window's contents; `folder` is the loaded results folder, `lib`
    /// holds the saved templates
    pub fn ui(&mut self, ui: &mut Ui, results: &Arc<Results>, folder: &str, lib: &mut Library) {
        // the first plot is in the default template's style
        if !self.started {
            self.started = true;
            if let Some(t) = lib.default_template() {
                self.style = t.apply(&self.style);
                self.template = t.name.clone();
                self.template_proteins = t.proteins;
            }
        }
        let results_id = Arc::as_ptr(results) as usize;
        let key = format!("{:?}", self.style.selection());
        let stale = !matches!(&self.data, Some((_, id, k)) if *id == results_id && *k == key);
        if stale {
            self.data = Some((Arc::new(circos::build(results, &self.style)), results_id, key));
            self.preview = None;
        }
        let data = self.data.as_ref().map(|(d, _, _)| d.clone()).unwrap_or_default();

        ui.horizontal_top(|ui| {
            egui::ScrollArea::vertical().id_salt("circos_controls").max_width(360.0).show(ui, |ui| {
                ui.set_width(340.0);
                // (a scroll area lays out like its parent: here, side by side)
                ui.vertical(|ui| self.controls(ui, &data, lib));
            });
            ui.separator();
            ui.vertical(|ui| {
                self.preview(ui, &data);
                self.export_bar(ui, &data, folder);
            });
        });
    }

    fn controls(&mut self, ui: &mut Ui, data: &CircosData, lib: &mut Library) {
        let s = &mut self.style;
        section(ui, "Crosslinks", true, |ui| {
            ui.horizontal(|ui| {
                ui.label("q-value ≤").on_hover_text("Both the CSM's spectrum q-value and its crosslink q-value.");
                ui.add(egui::DragValue::new(&mut s.max_q).range(0.0001..=1.0).speed(0.001).max_decimals(4));
                ui.label("CSMs ≥").on_hover_text("Only crosslinks seen in at least this many crosslink-spectrum matches.");
                ui.add(egui::DragValue::new(&mut s.min_csms).range(1..=1000));
            });
            ui.horizontal_wrapped(|ui| {
                ui.checkbox(&mut s.inter, "between proteins");
                ui.checkbox(&mut s.intra, "within a protein");
                ui.checkbox(&mut s.loops, "loop-links");
                ui.checkbox(&mut s.multivalent, "multivalent");
            });
            ui.checkbox(&mut s.all_places, "peptides of several proteins: in each").on_hover_text(
                "A peptide found in several proteins (or places in one) is linked in each of them; \
                 otherwise only in the first protein Sage lists for it.",
            );
            let csms: u32 = data.links.iter().map(|l| l.csms).sum();
            ui.label(RichText::new(format!("{} crosslinks · {csms} CSMs · {} proteins", data.links.len(), data.proteins.len())).weak());
        });

        section(ui, "Proteins", true, |ui| {
            ui.label(RichText::new("Clockwise from the start. Untick to leave a protein out; empty labels use the gene name.").weak());
            // the drawing order of every protein: the visible ones, then the hidden ones
            let mut order: Vec<String> = circos::shown_proteins(data, s).iter().map(|(p, _, _)| p.accession.clone()).collect();
            order.extend(data.proteins.iter().filter(|p| s.hidden.contains(&p.accession)).map(|p| p.accession.clone()));
            let mut swap: Option<(usize, usize)> = None;
            egui::Grid::new("circos_proteins").num_columns(5).spacing([6.0, 3.0]).show(ui, |ui| {
                for (i, accession) in order.iter().enumerate() {
                    let Some((index, p)) = data.proteins.iter().enumerate().find(|(_, p)| &p.accession == accession) else { continue };
                    let mut shown = !s.hidden.contains(accession);
                    if ui.checkbox(&mut shown, "").on_hover_text(accession.as_str()).changed() {
                        if shown {
                            s.hidden.retain(|a| a != accession);
                        } else {
                            s.hidden.push(accession.clone());
                        }
                    }
                    let label = s.labels.entry(accession.clone()).or_default();
                    // (a fixed size: a grid cell is only as wide as its column was)
                    ui.add_sized([110.0, ui.spacing().interact_size.y], egui::TextEdit::singleline(label).hint_text(&p.label));
                    let mut color = s.colors.get(accession).cloned().unwrap_or_else(|| s.palette.color(index).to_string());
                    if color_edit(ui, &mut color) {
                        s.colors.insert(accession.clone(), color);
                    }
                    ui.horizontal(|ui| {
                        if ui.add_enabled(i > 0, egui::Button::new("⏶").small()).clicked() {
                            swap = Some((i - 1, i));
                        }
                        if ui.add_enabled(i + 1 < order.len(), egui::Button::new("⏷").small()).clicked() {
                            swap = Some((i, i + 1));
                        }
                    });
                    ui.label(RichText::new(format!("{} links", p.links)).weak());
                    ui.end_row();
                }
            });
            s.labels.retain(|_, l| !l.trim().is_empty());
            if let Some((a, b)) = swap {
                order.swap(a, b);
                s.order = order;
            }
            ui.horizontal(|ui| {
                egui::ComboBox::from_id_salt("circos_palette").selected_text(s.palette.label()).show_ui(ui, |ui| {
                    for p in Palette::ALL {
                        ui.selectable_value(&mut s.palette, p, p.label());
                    }
                });
                if ui.button("Reset colours").clicked() {
                    s.colors.clear();
                }
                if ui.button("Reset order").clicked() {
                    s.order.clear();
                }
            });
        });

        section(ui, "Ring and labels", false, |ui| {
            number(ui, "start angle (°)", &mut s.start_angle, -360.0..=360.0, 1.0, "Where the first protein starts, clockwise from the top.");
            number(ui, "gap between proteins (°)", &mut s.gap, 0.0..=30.0, 0.1, "Space between neighbouring proteins.");
            number(ui, "ring radius", &mut s.radius, 0.3..=0.95, 0.005, "Outer radius of the protein ring, as a fraction of half the figure.");
            number(ui, "ring thickness", &mut s.thickness, 0.005..=0.2, 0.001, "As a fraction of half the figure.");
            ui.horizontal(|ui| {
                ui.label("protein labels");
                egui::ComboBox::from_id_salt("circos_labels").selected_text(s.label_mode.label()).show_ui(ui, |ui| {
                    for m in LabelMode::ALL {
                        ui.selectable_value(&mut s.label_mode, m, m.label());
                    }
                });
            });
            ui.horizontal(|ui| {
                ui.checkbox(&mut s.ticks, "ticks");
                ui.add_enabled(s.ticks, egui::Checkbox::new(&mut s.tick_labels, "numbered"));
                ui.label("every").on_hover_text("Residues between numbered ticks; 0 chooses from the longest protein.");
                ui.add_enabled(s.ticks, egui::DragValue::new(&mut s.tick_spacing).range(0..=100000).suffix(" res"));
            });
        });

        section(ui, "Tracks", false, |ui| {
            ui.label(RichText::new("Inside the ring, from the outside in.").weak());
            ui.checkbox(&mut s.coverage, "sequence coverage of the linear search").on_hover_text(
                "Residues covered by the linear search's target peptides at the q-value above.",
            );
            ui.horizontal(|ui| {
                ui.label("annotations").on_hover_text(
                    "A tab-separated file of regions, one per line: protein (accession, gene or entry name), \
                     start, end, label and an optional colour (#rrggbb). Lines starting with # are skipped.",
                );
                ui.add(egui::TextEdit::singleline(&mut s.annotations).desired_width(170.0).hint_text("domains.tsv"));
                if ui.button("Browse…").clicked() {
                    if let Some(path) = rfd::FileDialog::new().add_filter("tab-separated", &["tsv", "txt"]).pick_file() {
                        s.annotations = path.display().to_string();
                    }
                }
            });
            for note in &data.notes {
                ui.label(RichText::new(note).color(Color32::from_rgb(200, 120, 40)));
            }
            ui.horizontal(|ui| {
                ui.checkbox(&mut s.reactive, "reactive residues");
                ui.add_enabled(s.reactive, egui::TextEdit::singleline(&mut s.reactive_residues).desired_width(60.0).hint_text(&data.reactive))
                    .on_hover_text("Residues to mark; empty: the crosslinker's.");
            });
            ui.checkbox(&mut s.sites, "crosslinked residues (bar height: CSMs)");
            number(ui, "track width", &mut s.track_size, 0.005..=0.2, 0.001, "Each track's width, as a fraction of half the figure.");
        });

        section(ui, "Links", false, |ui| {
            ui.horizontal(|ui| {
                ui.label("colour by");
                egui::ComboBox::from_id_salt("circos_color_by").selected_text(s.color_by.label()).show_ui(ui, |ui| {
                    for c in ColorBy::ALL {
                        ui.selectable_value(&mut s.color_by, c, c.label());
                    }
                });
            });
            match s.color_by {
                ColorBy::Class => {
                    egui::Grid::new("circos_class_colors").num_columns(2).show(ui, |ui| {
                        for (name, color) in [
                            ("between proteins", &mut s.inter_color),
                            ("within a protein", &mut s.intra_color),
                            ("loop-links", &mut s.loop_color),
                            ("multivalent", &mut s.multi_color),
                        ] {
                            ui.label(name);
                            color_edit(ui, color);
                            ui.end_row();
                        }
                    });
                }
                ColorBy::Single => {
                    ui.horizontal(|ui| {
                        ui.label("colour");
                        color_edit(ui, &mut s.single_color);
                    });
                }
                ColorBy::Protein => {}
                _ => {
                    ui.horizontal(|ui| {
                        ui.label("low");
                        color_edit(ui, &mut s.low_color);
                        ui.label("high");
                        color_edit(ui, &mut s.high_color);
                    });
                }
            }
            number(ui, "line width (pt)", &mut s.link_width, 0.1..=10.0, 0.05, "Width of a link seen in one CSM.");
            ui.checkbox(&mut s.width_by_csms, "wider for more CSMs");
            number(ui, "opacity", &mut s.link_opacity, 0.05..=1.0, 0.01, "");
            number(ui, "bend (bezier radius)", &mut s.bezier_radius, 0.0..=1.0, 0.01, "How far the curves reach towards the centre: 0 through the centre, 1 straight.");
        });

        section(ui, "Highlights", false, |ui| highlights_ui(ui, s, &mut self.adding, data));

        section(ui, "Canvas", false, |ui| {
            number(ui, "size (inches)", &mut s.size, 2.0..=40.0, 0.1, "Width and height of the figure.");
            ui.horizontal(|ui| {
                ui.label("background");
                egui::ComboBox::from_id_salt("circos_background").selected_text(s.background.label()).show_ui(ui, |ui| {
                    for b in Background::ALL {
                        ui.selectable_value(&mut s.background, b, b.label());
                    }
                });
            });
            ui.horizontal(|ui| {
                ui.label("title");
                ui.add(egui::TextEdit::singleline(&mut s.title).desired_width(220.0));
            });
            ui.horizontal(|ui| {
                ui.checkbox(&mut s.legend, "legend");
                ui.add_enabled_ui(s.legend, |ui| {
                    egui::ComboBox::from_id_salt("circos_legend_corner").selected_text(s.legend_corner.label()).show_ui(ui, |ui| {
                        for c in LegendCorner::ALL {
                            ui.selectable_value(&mut s.legend_corner, c, c.label());
                        }
                    });
                })
                .response
                .on_hover_text("The legend's corner. The circle moves away from it, and is made smaller if it must, so that labels and legend never overlap.");
            });
        });

        section(ui, "Text", false, |ui| {
            let (fonts, filter) = (&mut self.fonts, &mut self.font_filter);
            ui.horizontal(|ui| {
                ui.label("font");
                let current = circos::first_font(&s.font);
                egui::ComboBox::from_id_salt("circos_font")
                    .width(230.0)
                    .height(380.0)
                    .close_behavior(egui::PopupCloseBehavior::CloseOnClickOutside)
                    .selected_text(current.clone())
                    .show_ui(ui, |ui| {
                        ui.add(TextEdit::singleline(filter).hint_text("search the installed fonts").desired_width(210.0));
                        let installed = fonts.get_or_insert_with(figure::font_families);
                        let needle = filter.trim().to_lowercase();
                        let mut picked = None;
                        for name in circos::GENERIC_FONTS.iter().map(|g| g.to_string()).chain(installed.iter().cloned()) {
                            if !needle.is_empty() && !name.to_lowercase().contains(&needle) {
                                continue;
                            }
                            if ui.selectable_label(name.eq_ignore_ascii_case(&current), name.as_str()).clicked() {
                                picked = Some(name);
                            }
                        }
                        if let Some(name) = picked {
                            s.font = circos::font_list(&name);
                            ui.memory_mut(|m| m.close_popup());
                        }
                    });
            });
            ui.horizontal(|ui| {
                ui.label("font list");
                ui.add(TextEdit::singleline(&mut s.font).desired_width(220.0)).on_hover_text(
                    "Fonts in order of preference, as in CSS: the first one installed is used. The SVG names \
                     them, so whoever opens it needs the font (or gets the next one); PDF and PNG carry it.",
                );
            });
            egui::Grid::new("circos_text").num_columns(4).spacing([8.0, 4.0]).show(ui, |ui| {
                ui.label("");
                ui.label(RichText::new("size (pt)").weak());
                ui.label("");
                ui.label("");
                ui.end_row();
                for (name, t) in [
                    ("title", &mut s.title_text),
                    ("protein labels", &mut s.label_text),
                    ("tick numbers", &mut s.tick_text),
                    ("legend", &mut s.legend_text),
                    ("annotation labels", &mut s.region_text),
                ] {
                    ui.label(name);
                    ui.add(egui::DragValue::new(&mut t.size).range(2.0..=96.0).speed(0.1).max_decimals(1));
                    ui.toggle_value(&mut t.bold, crate::theme::bold("B")).on_hover_text("bold");
                    ui.toggle_value(&mut t.italic, RichText::new("I").italics()).on_hover_text("italic");
                    ui.end_row();
                }
            });
            ui.horizontal(|ui| {
                ui.label("all text");
                if ui.button("smaller").on_hover_text("Every size 10% smaller").clicked() {
                    s.scale_text(1.0 / 1.1);
                }
                if ui.button("larger").on_hover_text("Every size 10% larger").clicked() {
                    s.scale_text(1.1);
                }
            });
        });

        section(ui, "Templates", true, |ui| self.templates(ui, data, lib));

        section(ui, "Built-in styles and style files", false, |ui| {
            ui.horizontal(|ui| {
                egui::ComboBox::from_id_salt("circos_preset").selected_text(self.preset.label()).show_ui(ui, |ui| {
                    for p in Preset::ALL {
                        ui.selectable_value(&mut self.preset, p, p.label());
                    }
                });
                if ui.button("Apply").on_hover_text("The preset's look; crosslinks, proteins, labels and text stay.").clicked() {
                    self.style = self.preset.apply(&self.style);
                }
            });
            ui.horizontal(|ui| {
                if ui
                    .button("Save style…")
                    .on_hover_text(
                        "The whole style as a JSON file: to share with others (Load style…, then save it as a \
                         template), or for `sage-gui --circos <results folder> <style.json>`.",
                    )
                    .clicked()
                {
                    if let Some(path) = rfd::FileDialog::new().set_file_name("circos_style.json").add_filter("JSON", &["json"]).save_file() {
                        self.message = Some(match std::fs::write(&path, self.style.to_json()) {
                            Ok(()) => (false, format!("Saved {}", path.display())),
                            Err(e) => (true, format!("Could not save {}: {e}", path.display())),
                        });
                    }
                }
                if ui.button("Load style…").clicked() {
                    if let Some(path) = rfd::FileDialog::new().add_filter("JSON", &["json"]).pick_file() {
                        match std::fs::read_to_string(&path).map_err(|e| e.to_string()).and_then(|t| CircosStyle::from_json(&t)) {
                            Ok(style) => {
                                self.style = style;
                                self.message = Some((false, format!("Loaded {}", path.display())));
                                // ready to be saved as a template of the same name
                                if self.template_name.trim().is_empty() {
                                    let stem = path.file_stem().map(|s| s.to_string_lossy().to_string()).unwrap_or_default();
                                    self.template_name = stem;
                                }
                            }
                            Err(e) => self.message = Some((true, format!("Could not load {}: {e}", path.display()))),
                        }
                    }
                }
                if ui.button("Reset").clicked() {
                    self.style = CircosStyle::default();
                }
            });
        });
    }

    /// The Templates section: the plot's style saved under a name, applied
    /// to other results
    fn templates(&mut self, ui: &mut Ui, data: &CircosData, lib: &mut Library) {
        ui.label(
            RichText::new(
                "Your styles, saved under a name on this computer. Apply one to other results to give their \
                 plot the same look.",
            )
            .weak(),
        );
        let picked = lib.template(&self.template).map(|(i, _)| i);
        ui.horizontal(|ui| {
            let text = match picked {
                Some(i) => template_label(lib, i),
                None if lib.circos_templates.is_empty() => "no templates yet".to_string(),
                None => "pick a template".to_string(),
            };
            egui::ComboBox::from_id_salt("circos_template").width(220.0).selected_text(text).show_ui(ui, |ui| {
                for i in 0..lib.circos_templates.len() {
                    let t = &lib.circos_templates[i];
                    if ui.selectable_label(picked == Some(i), template_label(lib, i)).on_hover_text(t.summary()).clicked() {
                        self.template = t.name.clone();
                        self.template_proteins = t.proteins;
                        self.confirm = None;
                        self.template_message = None;
                    }
                }
            });
            if ui
                .add_enabled(picked.is_some(), egui::Button::new("Apply"))
                .on_hover_text(
                    "The template's style for this plot. A template saved with its proteins also sets which \
                     are shown, their order, labels and colours, the annotation file and the title; one \
                     without leaves them as they are.",
                )
                .clicked()
            {
                if let Some(t) = picked.map(|i| &lib.circos_templates[i]) {
                    self.style = t.apply(&self.style);
                    self.template_message = Some((false, format!("Applied \"{}\".", t.name)));
                }
            }
        });
        if let Some(t) = picked.map(|i| &lib.circos_templates[i]) {
            ui.label(RichText::new(t.summary()).weak().small());
        }
        ui.horizontal(|ui| {
            let some = picked.is_some();
            if ui
                .add_enabled(some, egui::Button::new("Update"))
                .on_hover_text("Replace the template with the plot's style now (and its proteins, if ticked below).")
                .clicked()
            {
                self.confirm = Some(Confirm::Update);
            }
            if ui.add_enabled(some, egui::Button::new("Rename")).clicked() {
                self.confirm = Some(Confirm::Rename(self.template.clone()));
            }
            if ui.add_enabled(some, egui::Button::new("Delete")).clicked() {
                self.confirm = Some(Confirm::Delete);
            }
            let mut default = picked.is_some_and(|i| lib.is_default_template(i));
            if ui
                .add_enabled(some, egui::Checkbox::new(&mut default, "default"))
                .on_hover_text(
                    "New plots start from this template: this window when it opens, the plot written after \
                     each crosslink search, and `sage-gui --circos`. Marked ★ in the menu.",
                )
                .changed()
            {
                let done = if default { format!("New plots start from \"{}\".", self.template) } else { "New plots start from the built-in style.".into() };
                self.template_message = outcome(lib.set_default_template(if default { picked } else { None }), done);
            }
        });
        if let (Some(i), Some(confirm)) = (picked, self.confirm.clone()) {
            let name = lib.circos_templates[i].name.clone();
            ui.horizontal_wrapped(|ui| match confirm {
                Confirm::Update => {
                    let what = if self.template_proteins { "style and proteins" } else { "style" };
                    ui.label(format!("Replace \"{name}\" with the plot's {what} now?"));
                    if ui.button("Replace").clicked() {
                        let t = Template::new(&name, &self.style, self.template_proteins, data);
                        self.template_message = outcome(lib.put_template(t, Some(i)), format!("Updated \"{name}\"."));
                        self.confirm = None;
                    }
                    if ui.button("Cancel").clicked() {
                        self.confirm = None;
                    }
                }
                Confirm::Rename(mut new) => {
                    ui.label("New name");
                    let field = ui.add(TextEdit::singleline(&mut new).desired_width(150.0));
                    let enter = field.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
                    let rename = ui.button("Rename").clicked() || enter;
                    let cancel = ui.button("Cancel").clicked();
                    self.confirm = Some(Confirm::Rename(new.clone()));
                    if cancel {
                        self.confirm = None;
                    } else if rename {
                        let mut t = lib.circos_templates[i].clone();
                        t.name = new.trim().to_string();
                        match lib.put_template(t, Some(i)) {
                            Ok(()) => {
                                self.template_message = Some((false, format!("Renamed \"{name}\" to \"{}\".", new.trim())));
                                self.template = new.trim().to_string();
                                self.confirm = None;
                            }
                            Err(e) => self.template_message = Some((true, e)),
                        }
                    }
                }
                Confirm::Delete => {
                    ui.label(format!("Delete \"{name}\"? This cannot be undone."));
                    if ui.button("Delete").clicked() {
                        self.template_message = outcome(lib.delete_template(i), format!("Deleted \"{name}\"."));
                        self.template.clear();
                        self.confirm = None;
                    }
                    if ui.button("Cancel").clicked() {
                        self.confirm = None;
                    }
                }
            });
        }
        ui.horizontal(|ui| {
            ui.label("Save as");
            let field = ui.add(TextEdit::singleline(&mut self.template_name).hint_text("template name").desired_width(170.0));
            let enter = field.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)) && !self.template_name.trim().is_empty();
            if ui.button("Save").on_hover_text("A new template with the plot's style now.").clicked() || enter {
                let t = Template::new(&self.template_name, &self.style, self.template_proteins, data);
                let name = t.name.clone();
                match lib.put_template(t, None) {
                    Ok(()) => {
                        self.template_message = Some((false, format!("Saved \"{name}\".")));
                        self.template = name;
                        self.template_proteins = lib.template(&self.template).is_some_and(|(_, t)| t.proteins);
                        self.template_name.clear();
                        self.confirm = None;
                    }
                    Err(e) => self.template_message = Some((true, e)),
                }
            }
        });
        ui.checkbox(&mut self.template_proteins, "with these proteins").on_hover_text(
            "Save and Update also keep the proteins (matched by accession): which are shown, their order, \
             labels and colours (palette colours as they are now), the annotation file and the title. For \
             other results of the same proteins, e.g. a complex across experiments. Unticked, a template \
             leaves each plot's own proteins alone.",
        );
        if let Some(e) = &lib.error {
            show_message(ui, &Some((true, e.clone())));
        }
        show_message(ui, &self.template_message);
    }

    /// The figure as it will be exported, scaled to the space left
    fn preview(&mut self, ui: &mut Ui, data: &Arc<CircosData>) {
        let side = ui.available_width().min(ui.available_height() - 70.0).max(240.0);
        let pixels = ((side * ui.ctx().pixels_per_point()) / 40.0).round() as u32 * 40;
        let data_id = Arc::as_ptr(data) as usize;
        let current = matches!(&self.preview, Some((_, style, id, px)) if *style == self.style && *id == data_id && *px == pixels);
        if !current {
            let svg = circos::svg(data, &self.style);
            let dpi = ((pixels as f64 / self.style.size.clamp(2.0, 40.0)).round() as u32).max(10);
            match figure::rasterize(&svg, dpi) {
                Ok(image) => {
                    let size = [image.width() as usize, image.height() as usize];
                    let color = egui::ColorImage::from_rgba_unmultiplied(size, image.as_raw());
                    let texture = ui.ctx().load_texture("circos_preview", color, egui::TextureOptions::LINEAR);
                    self.preview = Some((texture, self.style.clone(), data_id, pixels));
                }
                Err(e) => {
                    self.message = Some((true, format!("Could not draw the plot: {e}")));
                    self.preview = None;
                }
            }
        }
        if let Some((texture, _, _, _)) = &self.preview {
            let frame = egui::Frame::default().fill(if self.style.background == Background::Black { Color32::BLACK } else { Color32::WHITE });
            frame.show(ui, |ui| {
                ui.add(egui::Image::new(texture).fit_to_exact_size(egui::vec2(side, side)));
            });
        }
    }

    fn export_bar(&mut self, ui: &mut Ui, data: &Arc<CircosData>, folder: &str) {
        ui.horizontal_wrapped(|ui| {
            ui.label(crate::theme::bold("Export"));
            egui::ComboBox::from_id_salt("circos_format").selected_text(self.format.label()).show_ui(ui, |ui| {
                for f in figure::Format::ALL {
                    ui.selectable_value(&mut self.format, f, f.label());
                }
            });
            if !self.format.is_vector() {
                ui.add(egui::DragValue::new(&mut self.dpi).range(72..=1200).suffix(" dpi"));
            }
            if ui.button("Save figure…").clicked() {
                let format = self.format;
                if let Some(path) = rfd::FileDialog::new()
                    .set_file_name(format!("crosslinks_circos.{}", format.extension()))
                    .add_filter(format.label(), &[format.extension()])
                    .save_file()
                {
                    let svg = circos::svg(data, &self.style);
                    let opaque = self.style.background != Background::Transparent;
                    self.message = Some(match figure::export(&svg, format, self.dpi, opaque, &path) {
                        Ok(()) => (false, format!("Saved {}", path.display())),
                        Err(e) => (true, format!("Could not save {}: {e}", path.display())),
                    });
                }
            }
            if ui
                .button("Circos files…")
                .on_hover_text(
                    "The plot as input files for Circos itself (circos.ca): karyotype, links, tracks and a \
                     circos.conf with this style. Run `circos -conf circos.conf` in the folder.",
                )
                .clicked()
            {
                if let Some(dir) = rfd::FileDialog::new().pick_folder() {
                    self.message = Some(match circos::write_circos_files(data, &self.style, &dir) {
                        Ok(files) => (false, format!("Wrote {} Circos files to {}", files.len(), dir.display())),
                        Err(e) => (true, format!("Could not write the Circos files: {e}")),
                    });
                }
            }
            if ui
                .button("Table…")
                .on_hover_text("Crosslinks between each pair of proteins, as a table for Circos' table viewer.")
                .clicked()
            {
                if let Some(path) = rfd::FileDialog::new().set_file_name("crosslink_table.txt").save_file() {
                    self.message = Some(match std::fs::write(&path, circos::crosslink_table(data, &self.style)) {
                        Ok(()) => (false, format!("Saved {}", path.display())),
                        Err(e) => (true, format!("Could not save {}: {e}", path.display())),
                    });
                }
            }
            if ui
                .add_enabled(!folder.is_empty(), egui::Button::new("All to results folder"))
                .on_hover_text(format!(
                    "SVG, PDF, PNG, this style, the Circos files and the table into {}",
                    Path::new(folder).join(circos::FOLDER).display()
                ))
                .clicked()
            {
                self.message = Some(match circos::write_figures(Path::new(folder), &self.style) {
                    Ok(files) => (false, format!("Wrote {} files to {}", files.len(), Path::new(folder).join(circos::FOLDER).display())),
                    Err(e) => (true, format!("Could not write the plots: {e}")),
                });
            }
        });
        if let Some((error, text)) = &self.message {
            let color = if *error { Color32::from_rgb(210, 60, 60) } else { ui.visuals().weak_text_color() };
            ui.label(RichText::new(text).color(color));
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use egui_kittest::{kittest::Queryable, Harness};

    /// The window on a real crosslink results folder (`SAGE_GUI_RESULTS`),
    /// with two saved templates, the first the default:
    /// `circos_window_*.png` in `SAGE_GUI_FIGURES`
    #[test]
    #[ignore]
    fn circos_window_preview() {
        let folder = std::env::var("SAGE_GUI_RESULTS").expect("set SAGE_GUI_RESULTS");
        let out = std::path::PathBuf::from(std::env::var("SAGE_GUI_FIGURES").expect("set SAGE_GUI_FIGURES"));
        std::fs::create_dir_all(&out).unwrap();
        let results = Arc::new(Results::load(Path::new(&folder)).unwrap());
        let mut lib = Library::default();
        let paper = CircosStyle {
            size: 6.0,
            palette: Palette::OkabeIto,
            font: "Times New Roman, serif".into(),
            label_text: circos::TypeStyle::new(10.0, true),
            gap: 3.0,
            ..CircosStyle::default()
        };
        lib.put_template(Template::new("Paper (Okabe-Ito, 6 in)", &paper, false, &CircosData::default()), None).unwrap();
        let data = circos::build(&results, &CircosStyle::default());
        let subunits = CircosStyle { title: "RNA polymerase II".into(), color_by: ColorBy::Protein, ..paper.clone() };
        lib.put_template(Template::new("Pol II subunits", &subunits, true, &data), None).unwrap();
        lib.set_default_template(Some(0)).unwrap();

        let variants = [
            ("light", false, None, false, false),
            ("dark", true, Some(Confirm::Update), false, false),
            ("text", false, None, true, false),
            ("highlights", false, None, false, true),
        ];
        let accession = |label: &str| data.proteins.iter().find(|p| p.label == label).map(|p| p.accession.clone()).unwrap_or_default();
        let strongest = data.links.iter().filter(|l| l.class == crate::viewer::mobiligram::Class::InterProtein).max_by_key(|l| l.csms).cloned();
        let marked = vec![
            Highlight::new(Pick::Between { a: accession("RPO21"), b: accession("RPB2") }, circos::HIGHLIGHT_COLORS[0]),
            Highlight::new(Pick::Protein { protein: accession("RPB3") }, circos::HIGHLIGHT_COLORS[1]),
            Highlight::new(Pick::Kind { kind: LinkKind::Loop }, circos::HIGHLIGHT_COLORS[3]),
        ];
        for (name, dark, confirm, text, highlights) in variants {
            let marked = marked.clone();
            let strongest = strongest.clone();
            let mut view = CircosView { open: true, ..Default::default() };
            let (results, folder) = (results.clone(), folder.clone());
            let mut lib = {
                let mut copy = Library::default();
                copy.circos_templates = lib.circos_templates.clone();
                copy.circos_default = lib.circos_default.clone();
                copy
            };
            let mut installed = false;
            let mut harness = Harness::builder()
                .with_size(egui::vec2(1240.0, 1500.0))
                .build(move |ctx| {
                    if !installed {
                        crate::theme::install(ctx, dark);
                        crate::theme::set_dark(ctx, dark);
                        installed = true;
                        return;
                    }
                    if let Some(c) = &confirm {
                        if view.started && view.confirm.is_none() {
                            view.confirm = Some(c.clone());
                            view.template_name = "Poster".into();
                        }
                    }
                    // highlights, the others dimmed; a crosslink on its way in
                    if highlights && view.started && view.style.highlights.is_empty() {
                        view.style.highlights = marked.clone();
                        view.style.dim_others = true;
                        view.adding.link = strongest.as_ref().map(|l| (l.a.clone(), l.b.clone()));
                    }
                    // italic protein labels, bold tick numbers, a larger legend
                    if text && view.started && !view.style.label_text.italic {
                        view.style.label_text.italic = true;
                        view.style.tick_text.bold = true;
                        view.style.legend_text.size = 10.0;
                    }
                    egui::Window::new("Crosslink plot (Circos)")
                        .default_size([1200.0, 1440.0])
                        .default_pos([10.0, 10.0])
                        .show(ctx, |ui| view.ui(ui, &results, &folder, &mut lib));
                });
            harness.run_steps(8);
            if highlights {
                harness.get_by_label("Proteins").click();
                harness.get_by_label("Highlights").click();
                harness.run_steps(4);
            }
            if text {
                // the Text section open (Proteins closed to make room), then the font menu
                harness.get_by_label("Proteins").click();
                harness.get_by_label("Text").click();
                harness.run_steps(4);
            }
            let path = out.join(format!("circos_window_{name}.png"));
            harness.render().expect("render").save(&path).unwrap();
            println!("{}", path.display());
            if text {
                // (the first menu of the controls, with Proteins closed)
                harness.get_all_by_role(egui::accesskit::Role::ComboBox).next().expect("font menu").click();
                harness.run_steps(4);
                let path = out.join("circos_window_fonts.png");
                harness.render().expect("render").save(&path).unwrap();
                println!("{}", path.display());
            }
        }
    }
}
