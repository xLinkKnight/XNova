//! The "Protein network" window of the Annotated spectra tab: the loaded
//! results' crosslinked proteins as an interaction network, what UniProt and
//! STRING know about them looked up on request, every choice of
//! [`NetworkStyle`] at hand, and a preview rendered from the very SVG that is
//! exported, in which proteins can be moved by dragging them.

use super::circos::{Background, Palette};
use super::circos_view::{color_edit, font_menu, number, section, show_message};
use super::figure;
use super::model::Results;
use super::network::{
    self, EdgeColor, EdgeWidth, Evidence, GroupBy, LabelPlace, Layout, LegendSide, NetworkData, NetworkStyle, NodeLabel, NodeSize, OwnGroups, Plan,
};
use super::protein_info::{self, Annotations};
use eframe::egui::{self, Color32, RichText, TextEdit, Ui};
use std::path::Path;
use std::sync::{mpsc, Arc};

/// What the lookup thread says
enum Lookup {
    Step(String),
    Done(Annotations, Vec<String>),
}

pub struct NetworkView {
    pub open: bool,
    style: NetworkStyle,
    /// The network's data, the results it came from and the selection it was built for
    data: Option<(Arc<NetworkData>, usize, String)>,
    /// What was looked up, for which results folder, and a count of changes
    ann: Arc<Annotations>,
    ann_folder: Option<String>,
    ann_version: u64,
    /// The user's groups, and the file they were read from
    own: (String, Arc<OwnGroups>),
    /// The plan, the key of its layout and the style it was made for
    plan: Option<(Arc<Plan>, String, NetworkStyle)>,
    /// The preview, and what it shows: style, plan and width in pixels
    preview: Option<(egui::TextureHandle, NetworkStyle, usize, u32)>,
    lookup: Option<mpsc::Receiver<Lookup>>,
    lookup_step: String,
    /// Look every protein up again, not only the new ones
    again: bool,
    lookup_message: Option<(bool, String)>,
    format: figure::Format,
    dpi: u32,
    message: Option<(bool, String)>,
    fonts: Option<Vec<String>>,
    font_filter: String,
    /// The protein being dragged in the preview
    dragging: Option<String>,
    /// Where the preview was last drawn
    preview_rect: Option<egui::Rect>,
}

impl Default for NetworkView {
    fn default() -> Self {
        NetworkView {
            open: false,
            style: NetworkStyle::default(),
            data: None,
            ann: Arc::new(Annotations::default()),
            ann_folder: None,
            ann_version: 0,
            own: (String::new(), Arc::new(OwnGroups::default())),
            plan: None,
            preview: None,
            lookup: None,
            lookup_step: String::new(),
            again: false,
            lookup_message: None,
            format: figure::Format::Svg,
            dpi: 300,
            message: None,
            fonts: None,
            font_filter: String::new(),
            dragging: None,
            preview_rect: None,
        }
    }
}

/// A combo box of an enum's values
fn choice<T: Copy + PartialEq>(ui: &mut Ui, id: &str, label: &str, value: &mut T, all: &[T], name: impl Fn(T) -> &'static str, width: f32) {
    ui.horizontal(|ui| {
        ui.label(label);
        egui::ComboBox::from_id_salt(id).width(width).selected_text(name(*value)).show_ui(ui, |ui| {
            for &v in all {
                ui.selectable_value(value, v, name(v));
            }
        });
    });
}

/// What is known about a protein, for its tooltip
fn about(ui: &mut Ui, node: &network::Node, info: Option<&protein_info::ProteinInfo>, label: &str, group: Option<&str>) {
    ui.set_max_width(420.0);
    ui.label(crate::theme::bold(label));
    let name = info.map(|p| p.name.as_str()).filter(|n| !n.is_empty()).unwrap_or(&node.name);
    if !name.is_empty() {
        ui.label(name);
    }
    ui.label(RichText::new(&node.accession).weak());
    if let Some(p) = info {
        if !p.organism.is_empty() {
            ui.label(RichText::new(&p.organism).italics().weak());
        }
    }
    ui.label(format!(
        "{} crosslink{} to {} protein{}, {} within it; {} CSMs",
        node.crosslinks,
        if node.crosslinks == 1 { "" } else { "s" },
        node.partners,
        if node.partners == 1 { "" } else { "s" },
        node.self_links,
        node.csms + node.self_csms
    ));
    if let Some(g) = group {
        ui.label(format!("group: {g}"));
    }
    let Some(p) = info else {
        ui.label(RichText::new("Not looked up (or not in UniProt).").weak());
        return;
    };
    let list = |terms: &[protein_info::Term], most: usize| {
        let mut text = terms.iter().take(most).map(|t| t.name.as_str()).collect::<Vec<_>>().join("; ");
        if terms.len() > most {
            text.push_str(&format!("; and {} more", terms.len() - most));
        }
        text
    };
    let mut row = |title: &str, text: String| {
        if !text.is_empty() {
            ui.horizontal_wrapped(|ui| {
                ui.label(crate::theme::bold(title));
                ui.label(text);
            });
        }
    };
    row("location", p.locations.join("; "));
    row("complexes", list(&p.complexes, 6));
    row("pathways", list(&p.pathways, 5));
    row("GO process", list(&p.go_process, 5));
    row("GO function", list(&p.go_function, 4));
    if !p.function.is_empty() {
        let mut f: String = p.function.chars().take(420).collect();
        if p.function.chars().count() > 420 {
            f.push('…');
        }
        row("function", f);
    }
}

impl NetworkView {
    /// The window's contents; `folder` is the loaded results folder
    pub fn ui(&mut self, ui: &mut Ui, results: &Arc<Results>, folder: &str) {
        // what was looked up for these results (kept in their folder)
        if self.ann_folder.as_deref() != Some(folder) {
            self.ann = Arc::new(if folder.is_empty() { None } else { Annotations::load(Path::new(folder)) }.unwrap_or_default());
            self.ann_folder = Some(folder.to_string());
            self.ann_version += 1;
            self.lookup_message = None;
        }
        self.poll_lookup(ui.ctx(), folder);

        let results_id = Arc::as_ptr(results) as usize;
        let key = format!("{:?}", self.style.selection());
        if !matches!(&self.data, Some((_, id, k)) if *id == results_id && *k == key) {
            self.data = Some((Arc::new(network::build(results, &self.style)), results_id, key));
            self.plan = None;
        }
        let data = self.data.as_ref().map(|(d, _, _)| d.clone()).unwrap_or_default();
        if self.own.0 != self.style.own_groups {
            self.own = (self.style.own_groups.clone(), Arc::new(OwnGroups::read(&self.style.own_groups)));
        }
        let plan = self.plan_for(&data);

        ui.horizontal_top(|ui| {
            egui::ScrollArea::vertical().id_salt("network_controls").max_width(360.0).show(ui, |ui| {
                ui.set_width(340.0);
                ui.vertical(|ui| self.controls(ui, &data, &plan, folder));
            });
            ui.separator();
            ui.vertical(|ui| {
                self.preview(ui, &data, &plan);
                self.export_bar(ui, &data, &plan, folder);
            });
        });
    }

    /// The plan for the current style: the layout kept while only colours
    /// and labels of groups change
    fn plan_for(&mut self, data: &Arc<NetworkData>) -> Arc<Plan> {
        let layout_key = format!("{}|{}|{}", Arc::as_ptr(data) as usize, self.ann_version, self.style.plan_key());
        match &self.plan {
            Some((plan, key, style)) if *key == layout_key && *style == self.style => plan.clone(),
            previous => {
                let layout = previous.as_ref().filter(|(_, key, _)| *key == layout_key).map(|(p, _, _)| (p.layout.clone(), p.aspect));
                let plan = Arc::new(network::plan_with(data, &self.ann, &self.own.1, &self.style, layout));
                self.plan = Some((plan.clone(), layout_key, self.style.clone()));
                plan
            }
        }
    }

    fn poll_lookup(&mut self, ctx: &egui::Context, folder: &str) {
        let Some(rx) = &self.lookup else { return };
        loop {
            match rx.try_recv() {
                Ok(Lookup::Step(step)) => self.lookup_step = step,
                Ok(Lookup::Done(ann, problems)) => {
                    let kept = if folder.is_empty() {
                        Err("no results folder".to_string())
                    } else {
                        ann.save(Path::new(folder))
                    };
                    let found = ann.proteins.len();
                    self.ann = Arc::new(ann);
                    self.ann_version += 1;
                    self.lookup = None;
                    self.lookup_step.clear();
                    let mut text = format!("{found} proteins known");
                    match kept {
                        Ok(path) => text.push_str(&format!("; kept in {}", path.display())),
                        Err(e) => text.push_str(&format!("; could not be kept ({e}), used until the window closes")),
                    }
                    self.lookup_message = Some(if problems.is_empty() { (false, text) } else { (true, format!("{}. {text}", problems.join(". "))) });
                    return;
                }
                Err(mpsc::TryRecvError::Empty) => {
                    ctx.request_repaint_after(std::time::Duration::from_millis(200));
                    return;
                }
                Err(mpsc::TryRecvError::Disconnected) => {
                    self.lookup = None;
                    self.lookup_message = Some((true, "The lookup stopped unexpectedly".into()));
                    return;
                }
            }
        }
    }

    fn start_lookup(&mut self, ctx: &egui::Context, data: &NetworkData) {
        let (tx, rx) = mpsc::channel();
        let accessions = data.accessions();
        let have = (*self.ann).clone();
        let again = self.again;
        let ctx = ctx.clone();
        let started = std::thread::Builder::new().name("protein lookup".into()).spawn(move || {
            let progress = |step: String| {
                let _ = tx.send(Lookup::Step(step));
                ctx.request_repaint();
            };
            let (ann, problems) = protein_info::look_up(&accessions, &have, again, &progress);
            let _ = tx.send(Lookup::Done(ann, problems));
            ctx.request_repaint();
        });
        match started {
            Ok(_) => {
                self.lookup = Some(rx);
                self.lookup_step = "starting…".into();
                self.lookup_message = None;
            }
            Err(e) => self.lookup_message = Some((true, format!("Could not start the lookup: {e}"))),
        }
    }

    fn controls(&mut self, ui: &mut Ui, data: &Arc<NetworkData>, plan: &Plan, folder: &str) {
        // what is known about the proteins
        section(ui, "Protein information", true, |ui| {
            let ann = self.ann.clone();
            let total = data.nodes.len();
            if ann.is_empty() {
                ui.label(
                    RichText::new(
                        "Not looked up yet. Look the proteins up to group them by complex, pathway, GO term or \
                         location, and to tell known interactions from new ones.",
                    )
                    .weak(),
                );
            } else {
                let found = data.nodes.iter().filter(|n| ann.proteins.contains_key(&n.accession)).count();
                let mut text = format!("{found} of {total} proteins found in UniProt");
                if !ann.uniprot_release.is_empty() {
                    text.push_str(&format!(" (release {})", ann.uniprot_release));
                }
                if !ann.string_version.is_empty() {
                    text.push_str(&format!("; their interactions from STRING {}", ann.string_version));
                }
                if !ann.looked_up.is_empty() {
                    text.push_str(&format!("; looked up {}", ann.looked_up));
                }
                ui.label(text);
                let missing: Vec<&String> = data.nodes.iter().map(|n| &n.accession).filter(|a| ann.not_found.contains(*a)).collect();
                if !missing.is_empty() {
                    ui.label(RichText::new(format!("{} not found (hover to list)", missing.len())).weak()).on_hover_text(
                        missing.iter().map(|a| a.as_str()).collect::<Vec<_>>().join("\n")
                            + "\n\nUniProt accessions only (sp|P04050|RPB1_YEAST, P04050); for other proteins, \
                               group them with a file of your own (Groups).",
                    );
                }
                let new = ann.missing(&data.accessions()).len();
                if new > 0 {
                    ui.label(RichText::new(format!("{new} not looked up yet (new in these results or settings)")).weak());
                }
            }
            ui.horizontal(|ui| {
                let running = self.lookup.is_some();
                if ui
                    .add_enabled(!running && total > 0, egui::Button::new("Look up in UniProt and STRING"))
                    .on_hover_text(
                        "UniProt: name, gene, organism, function, subcellular location, GO terms, complexes \
                         (ComplexPortal), pathways (Reactome) and keywords. STRING: the known and predicted \
                         interactions between these proteins, with their evidence.",
                    )
                    .clicked()
                {
                    self.start_lookup(ui.ctx(), data);
                }
                ui.checkbox(&mut self.again, "again").on_hover_text("Every protein again, not only those not looked up yet");
            });
            if self.lookup.is_some() {
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label(RichText::new(&self.lookup_step).weak());
                });
            }
            show_message(ui, &self.lookup_message);
            ui.label(
                RichText::new(format!(
                    "Sends the proteins' accessions, and nothing else, to rest.uniprot.org and string-db.org. \
                     What comes back is kept in {} in the results folder, for the network to be drawn again offline.",
                    Path::new(network::FOLDER).join(protein_info::FILE).display()
                ))
                .small()
                .weak(),
            );
        });

        let s = &mut self.style;
        section(ui, "Crosslinks", false, |ui| {
            ui.horizontal(|ui| {
                ui.label("q-value at most").on_hover_text("CSM and unique-crosslink q-value, as the crosslink plot");
                ui.add(egui::DragValue::new(&mut s.max_q).range(0.0..=1.0).speed(0.001).max_decimals(3));
                ui.label("CSMs at least");
                ui.add(egui::DragValue::new(&mut s.min_csms).range(1..=1000));
            });
            ui.checkbox(&mut s.intra, "crosslinks within a protein (a loop on its node)");
            ui.checkbox(&mut s.loops, "loop-links (a loop on its node)");
            ui.checkbox(&mut s.multivalent, "multivalent species (each of their links)");
            ui.checkbox(&mut s.all_places, "shared peptides in every protein")
                .on_hover_text("A peptide found in several proteins links each of them; otherwise only the first Sage listed");
        });

        section(ui, "Proteins", false, |ui| {
            choice(ui, "net_label_by", "label", &mut s.label_by, &NodeLabel::ALL, NodeLabel::label, 130.0);
            choice(ui, "net_label_place", "placed", &mut s.label_place, &LabelPlace::ALL, LabelPlace::label, 130.0);
            choice(ui, "net_size_by", "size by", &mut s.size_by, &NodeSize::ALL, NodeSize::label, 130.0);
            ui.horizontal(|ui| {
                ui.label("diameter (pt)");
                ui.add(egui::DragValue::new(&mut s.node_size).range(1.0..=200.0).speed(0.2).max_decimals(1));
                ui.label("to");
                ui.add(egui::DragValue::new(&mut s.node_size_max).range(1.0..=200.0).speed(0.2).max_decimals(1));
            });
            ui.horizontal(|ui| {
                ui.label("outline");
                color_edit(ui, &mut s.node_outline);
            });
            ui.checkbox(&mut s.lone, "proteins crosslinked only within themselves").on_hover_text("As nodes without edges");
            egui::CollapsingHeader::new(format!("Each protein ({})", data.nodes.len())).id_salt("net_each").show(ui, |ui| {
                ui.label(RichText::new("Untick to leave a protein out; type a label of your own.").weak());
                egui::ScrollArea::vertical().id_salt("net_protein_list").max_height(300.0).show(ui, |ui| {
                    egui::Grid::new("net_proteins").num_columns(3).spacing([6.0, 3.0]).show(ui, |ui| {
                        for n in &data.nodes {
                            let mut shown = !s.hidden.contains(&n.accession);
                            if ui.checkbox(&mut shown, "").changed() {
                                if shown {
                                    s.hidden.retain(|h| h != &n.accession);
                                } else {
                                    s.hidden.push(n.accession.clone());
                                }
                            }
                            let info = self.ann.proteins.get(&n.accession);
                            let hint = network::default_label(n, info, s.label_by);
                            let mut text = s.labels.get(&n.accession).cloned().unwrap_or_default();
                            if ui.add_sized([110.0, ui.spacing().interact_size.y], TextEdit::singleline(&mut text).hint_text(&hint)).changed() {
                                if text.trim().is_empty() {
                                    s.labels.remove(&n.accession);
                                } else {
                                    s.labels.insert(n.accession.clone(), text);
                                }
                            }
                            ui.label(RichText::new(format!("{} · {} xl", n.accession, n.crosslinks + n.self_links)).weak())
                                .on_hover_ui(|ui| about(ui, n, info, &hint, None));
                            ui.end_row();
                        }
                    });
                });
            });
        });

        section(ui, "Groups", true, |ui| {
            ui.horizontal(|ui| {
                ui.label("group by");
                egui::ComboBox::from_id_salt("net_group_by").width(220.0).selected_text(s.group_by.label()).show_ui(ui, |ui| {
                    for g in GroupBy::ALL {
                        ui.selectable_value(&mut s.group_by, g, g.label());
                    }
                });
            });
            if s.group_by.looked_up() && self.ann.proteins.is_empty() {
                ui.label(RichText::new("Needs the proteins looked up (Protein information).").color(Color32::from_rgb(200, 120, 40)));
            }
            if s.group_by == GroupBy::Own {
                ui.horizontal(|ui| {
                    ui.label("file");
                    ui.add(TextEdit::singleline(&mut s.own_groups).desired_width(190.0)).on_hover_text(
                        "A tab- or comma-separated file: a protein (accession, gene or entry name) and its group per \
                         line, optionally a colour (#rrggbb); lines starting with # are skipped.",
                    );
                    if ui.button("Choose…").clicked() {
                        if let Some(path) = rfd::FileDialog::new().add_filter("table", &["tsv", "txt", "csv"]).pick_file() {
                            s.own_groups = path.display().to_string();
                        }
                    }
                });
                for note in &self.own.1.notes {
                    ui.label(RichText::new(note).color(Color32::from_rgb(200, 120, 40)));
                }
            }
            if s.group_by == GroupBy::Module {
                if let Some(q) = plan.modularity {
                    ui.label(RichText::new(format!("modularity {q:.2} (from −0.5 to 1; above 0.3, clearly divided)")).weak());
                }
            }
            if s.group_by != GroupBy::None {
                ui.horizontal(|ui| {
                    ui.label("at least");
                    ui.add(egui::DragValue::new(&mut s.min_group).range(1..=100));
                    ui.label("proteins, at most");
                    ui.add(egui::DragValue::new(&mut s.max_groups).range(1..=60));
                    ui.label("groups");
                });
                ui.label(
                    RichText::new(
                        "Each protein joins the group shared by the most proteins, then the next, and so on; \
                         leave a group out (e.g. too general) for its proteins to join their next.",
                    )
                    .weak(),
                );
                if !plan.groups.is_empty() {
                    egui::Grid::new("net_groups").num_columns(4).spacing([6.0, 3.0]).show(ui, |ui| {
                        for g in &plan.groups {
                            let mut color = g.color.clone();
                            if color_edit(ui, &mut color) {
                                s.group_colors.insert(g.key.clone(), color);
                            }
                            let mut text = s.group_labels.get(&g.key).cloned().unwrap_or_default();
                            if ui
                                .add_sized([170.0, ui.spacing().interact_size.y], TextEdit::singleline(&mut text).hint_text(&g.name))
                                .on_hover_text(format!("{} ({}): the label (empty: its name)", g.name, g.key))
                                .changed()
                            {
                                if text.trim().is_empty() {
                                    s.group_labels.remove(&g.key);
                                } else {
                                    s.group_labels.insert(g.key.clone(), text);
                                }
                            }
                            ui.label(RichText::new(g.members.to_string()).weak()).on_hover_text("proteins");
                            if ui.add(egui::Button::new("leave out").small()).on_hover_text("Its proteins join their next group").clicked() {
                                s.preferred.retain(|p| p != &g.key);
                                s.excluded.push(g.key.clone());
                            }
                            ui.end_row();
                        }
                    });
                }
                if !s.excluded.is_empty() {
                    ui.label(RichText::new("left out").weak());
                    let mut back = None;
                    for (i, key) in s.excluded.iter().enumerate() {
                        ui.horizontal(|ui| {
                            if ui.add(egui::Button::new("use again").small()).clicked() {
                                back = Some(i);
                            }
                            ui.label(RichText::new(key).weak());
                        });
                    }
                    if let Some(i) = back {
                        s.excluded.remove(i);
                    }
                }
                if !s.preferred.is_empty() {
                    ui.horizontal_wrapped(|ui| {
                        ui.label(RichText::new(format!("preferred first: {}", s.preferred.join(", "))).weak());
                        if ui.add(egui::Button::new("clear").small()).clicked() {
                            s.preferred.clear();
                        }
                    });
                }
                if !plan.others.is_empty() {
                    egui::CollapsingHeader::new(format!("Other terms ({})", plan.others.len())).id_salt("net_others").show(ui, |ui| {
                        ui.label(
                            RichText::new(
                                "Proteins having each, among those drawn: their proteins are in a group taken first, \
                                 or too few. Prefer one to group by it first (e.g. the pathway to name a group by, \
                                 of several the same proteins share).",
                            )
                            .weak(),
                        );
                        egui::ScrollArea::vertical().id_salt("net_others_list").max_height(220.0).show(ui, |ui| {
                            for (key, name, count) in plan.others.iter().take(300) {
                                ui.horizontal(|ui| {
                                    if ui.add(egui::Button::new("prefer").small()).on_hover_text("Group by it first").clicked() {
                                        s.preferred.retain(|p| p != key);
                                        s.preferred.insert(0, key.clone());
                                    }
                                    ui.label(RichText::new(format!("{count} · {name}")).weak()).on_hover_text(key.as_str());
                                });
                            }
                        });
                    });
                }
                ui.horizontal(|ui| {
                    ui.label("palette");
                    egui::ComboBox::from_id_salt("net_palette").width(190.0).selected_text(s.palette.label()).show_ui(ui, |ui| {
                        for p in Palette::ALL {
                            ui.selectable_value(&mut s.palette, p, p.label());
                        }
                    });
                    if ui.add_enabled(!s.group_colors.is_empty(), egui::Button::new("reset").small()).on_hover_text("The palette's colours for every group").clicked() {
                        s.group_colors.clear();
                    }
                });
            }
            ui.horizontal(|ui| {
                ui.label(if s.group_by == GroupBy::None { "colour" } else { "in no group" });
                color_edit(ui, &mut s.ungrouped_color);
            });
            ui.checkbox(&mut s.group_names, "group names in the figure, next to their proteins");
        });

        section(ui, "Interactions (edges)", false, |ui| {
            choice(ui, "net_color_by", "colour", &mut s.color_by, &EdgeColor::ALL, EdgeColor::label, 190.0);
            match s.color_by {
                EdgeColor::Known => {
                    ui.horizontal(|ui| {
                        ui.label("known");
                        color_edit(ui, &mut s.known_color);
                        ui.label("new");
                        color_edit(ui, &mut s.new_color);
                    });
                }
                EdgeColor::Single => {
                    ui.horizontal(|ui| {
                        ui.label("colour");
                        color_edit(ui, &mut s.edge_color);
                    });
                }
                _ => {
                    ui.horizontal(|ui| {
                        ui.label("sharing none");
                        color_edit(ui, &mut s.other_color);
                    });
                }
            }
            choice(ui, "net_evidence", "known by", &mut s.evidence, &Evidence::ALL, Evidence::label, 190.0);
            ui.horizontal(|ui| {
                ui.label("STRING score at least");
                ui.add(egui::Slider::new(&mut s.min_score, 0.15..=0.99).max_decimals(2)).on_hover_text(
                    "STRING's confidence: 0.15 low, 0.4 medium, 0.7 high, 0.9 highest. Experiments: physical \
                     interaction data (BioGRID, IntAct and others); curated databases: complexes and pathways.",
                );
            });
            choice(ui, "net_width_by", "width by", &mut s.width_by, &EdgeWidth::ALL, EdgeWidth::label, 190.0);
            ui.horizontal(|ui| {
                ui.label("width (pt)");
                ui.add(egui::DragValue::new(&mut s.edge_width).range(0.05..=30.0).speed(0.05).max_decimals(2));
                ui.label("to");
                ui.add(egui::DragValue::new(&mut s.edge_width_max).range(0.05..=30.0).speed(0.05).max_decimals(2));
            });
            number(ui, "opacity", &mut s.edge_opacity, 0.05..=1.0, 0.01, "Of the edges");
            ui.checkbox(&mut s.numbers, "the number of crosslinks on each edge");
            ui.horizontal(|ui| {
                ui.label("loops (within a protein)");
                color_edit(ui, &mut s.self_color);
            });
        });

        section(ui, "Layout", false, |ui| {
            choice(ui, "net_layout", "layout", &mut s.layout, &Layout::ALL, Layout::label, 160.0);
            number(ui, "spacing", &mut s.spacing, 0.2..=5.0, 0.01, "Distance between the proteins, times the default");
            if s.layout != Layout::Circle {
                number(ui, "groups gather", &mut s.group_pull, 0.0..=5.0, 0.01, "How strongly the proteins of a group are drawn together (0: not at all)");
            }
            ui.horizontal(|ui| {
                if ui.button("Another layout").on_hover_text(format!("A layout from another start (now {})", s.seed)).clicked() {
                    s.seed = s.seed.wrapping_add(1);
                }
                if ui
                    .add_enabled(!s.positions.is_empty(), egui::Button::new(format!("Put moved proteins back ({})", s.positions.len())))
                    .clicked()
                {
                    s.positions.clear();
                }
            });
            ui.label(RichText::new("Drag a protein in the preview to move it; hover over one for what is known about it.").weak());
        });

        section(ui, "Canvas", false, |ui| {
            ui.horizontal(|ui| {
                ui.label("width");
                ui.add(egui::DragValue::new(&mut s.width).range(1.0..=40.0).speed(0.05).max_decimals(2).suffix(" in"));
                ui.label("height");
                ui.add(egui::DragValue::new(&mut s.height).range(1.0..=40.0).speed(0.05).max_decimals(2).suffix(" in"));
            });
            ui.horizontal(|ui| {
                ui.label("background");
                egui::ComboBox::from_id_salt("net_background").selected_text(s.background.label()).show_ui(ui, |ui| {
                    for b in Background::ALL {
                        ui.selectable_value(&mut s.background, b, b.label());
                    }
                });
            });
            ui.horizontal(|ui| {
                ui.label("title");
                ui.add(TextEdit::singleline(&mut s.title).desired_width(230.0));
            });
            ui.horizontal(|ui| {
                ui.checkbox(&mut s.legend, "legend");
                egui::ComboBox::from_id_salt("net_legend_side").width(80.0).selected_text(s.legend_side.label()).show_ui(ui, |ui| {
                    for side in LegendSide::ALL {
                        ui.selectable_value(&mut s.legend_side, side, side.label());
                    }
                });
                ui.checkbox(&mut s.summary, "with the numbers");
            });
        });

        section(ui, "Text", false, |ui| {
            font_menu(ui, "net_font", &mut s.font, &mut self.fonts, &mut self.font_filter);
            egui::Grid::new("net_text").num_columns(4).spacing([8.0, 4.0]).show(ui, |ui| {
                ui.label("");
                ui.label(RichText::new("size (pt)").weak());
                ui.label("");
                ui.label("");
                ui.end_row();
                for (name, t) in [
                    ("title", &mut s.title_text),
                    ("protein labels", &mut s.label_text),
                    ("group names", &mut s.group_text),
                    ("legend", &mut s.legend_text),
                    ("edge numbers", &mut s.number_text),
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

        section(ui, "Style files", false, |ui| {
            ui.horizontal(|ui| {
                if ui
                    .button("Save style…")
                    .on_hover_text("The whole style as a JSON file (with proteins moved by hand), to load again or share")
                    .clicked()
                {
                    if let Some(path) = rfd::FileDialog::new().set_file_name("network_style.json").save_file() {
                        self.message = Some(match std::fs::write(&path, s.to_json()) {
                            Ok(()) => (false, format!("Saved {}", path.display())),
                            Err(e) => (true, format!("Could not save {}: {e}", path.display())),
                        });
                    }
                }
                if ui.button("Load style…").clicked() {
                    if let Some(path) = rfd::FileDialog::new().add_filter("JSON", &["json"]).pick_file() {
                        match std::fs::read_to_string(&path).map_err(|e| e.to_string()).and_then(|t| NetworkStyle::from_json(&t)) {
                            Ok(style) => {
                                *s = style;
                                self.message = Some((false, format!("Loaded {}", path.display())));
                            }
                            Err(e) => self.message = Some((true, format!("{}: {e}", path.display()))),
                        }
                    }
                }
                if ui.button("Default style").on_hover_text("Every setting as it was first").clicked() {
                    *s = NetworkStyle::default();
                }
            });
            if !folder.is_empty() {
                let saved = Path::new(folder).join(network::FOLDER).join("network_style.json");
                if saved.is_file() && ui.button("Style of the results folder").on_hover_text(saved.display().to_string()).clicked() {
                    match std::fs::read_to_string(&saved).map_err(|e| e.to_string()).and_then(|t| NetworkStyle::from_json(&t)) {
                        Ok(style) => *s = style,
                        Err(e) => self.message = Some((true, format!("{}: {e}", saved.display()))),
                    }
                }
            }
        });
    }

    /// The figure as it will be exported, scaled to the space left; proteins
    /// are dragged to move them, and tell what is known about them on hover
    fn preview(&mut self, ui: &mut Ui, data: &Arc<NetworkData>, plan: &Arc<Plan>) {
        let (w_in, h_in) = (self.style.width.clamp(1.0, 40.0), self.style.height.clamp(1.0, 40.0));
        let room = egui::vec2(ui.available_width().max(240.0), (ui.available_height() - 70.0).max(240.0));
        let scale = (room.x / w_in as f32).min(room.y / h_in as f32);
        let size = egui::vec2(w_in as f32 * scale, h_in as f32 * scale);
        let pixels = ((size.x * ui.ctx().pixels_per_point()) / 40.0).round() as u32 * 40;
        let plan_id = Arc::as_ptr(plan) as usize;
        let current = matches!(&self.preview, Some((_, style, id, px)) if *style == self.style && *id == plan_id && *px == pixels);
        if !current {
            let svg = network::svg(data, &self.ann, plan, &self.style);
            let dpi = ((pixels as f64 / w_in).round() as u32).max(10);
            match figure::rasterize(&svg, dpi) {
                Ok(image) => {
                    let dims = [image.width() as usize, image.height() as usize];
                    let color = egui::ColorImage::from_rgba_unmultiplied(dims, image.as_raw());
                    let texture = ui.ctx().load_texture("network_preview", color, egui::TextureOptions::LINEAR);
                    self.preview = Some((texture, self.style.clone(), plan_id, pixels));
                }
                Err(e) => {
                    self.message = Some((true, format!("Could not draw the network: {e}")));
                    self.preview = None;
                }
            }
        }
        let Some((texture, _, _, _)) = &self.preview else { return };
        let frame = egui::Frame::default().fill(if self.style.background == Background::Black { Color32::BLACK } else { Color32::WHITE });
        let response = frame
            .show(ui, |ui| ui.add(egui::Image::new(texture).fit_to_exact_size(size).sense(egui::Sense::click_and_drag())))
            .inner;

        // the figure's units under the pointer, and the protein there
        let f = network::frame(data, &self.ann, plan, &self.style);
        let rect = response.rect;
        self.preview_rect = Some(rect);
        let to_units = |p: egui::Pos2| (((p.x - rect.min.x) / rect.width()) as f64 * f.width, ((p.y - rect.min.y) / rect.height()) as f64 * f.height);
        let units_per_point = f.width / rect.width() as f64;
        let hit = |p: egui::Pos2| -> Option<usize> {
            let (x, y) = to_units(p);
            f.nodes
                .iter()
                .enumerate()
                .map(|(k, n)| (k, (n.0 - x).hypot(n.1 - y) - n.2))
                .filter(|(_, d)| *d <= 4.0 * units_per_point)
                .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
                .map(|(k, _)| k)
        };
        let pointer = response.interact_pointer_pos().or(response.hover_pos());
        if response.drag_started() {
            self.dragging = pointer.and_then(hit).map(|k| data.nodes[plan.nodes[k]].accession.clone());
        }
        if let (Some(accession), Some(p)) = (&self.dragging, pointer) {
            if response.dragged() {
                // in the layout's box, anywhere in the area
                let (x, y) = to_units(p);
                let (x, y) = (x.clamp(f.area[0], f.area[2]), y.clamp(f.area[1], f.area[3]));
                let b = f.layout_box;
                let nx = (x - b[0]) / (b[2] - b[0]).max(1e-9);
                let ny = (y - b[1]) / (b[3] - b[1]).max(1e-9);
                self.style.positions.insert(accession.clone(), [nx, ny]);
            }
        }
        if response.drag_stopped() {
            self.dragging = None;
        }
        if self.dragging.is_some() {
            ui.ctx().set_cursor_icon(egui::CursorIcon::Grabbing);
        } else if let Some(k) = response.hover_pos().and_then(hit) {
            ui.ctx().set_cursor_icon(egui::CursorIcon::Grab);
            let node = &data.nodes[plan.nodes[k]];
            let info = self.ann.proteins.get(&node.accession);
            let group = plan.group_of[k].map(|g| plan.groups[g].label.as_str());
            response.on_hover_ui_at_pointer(|ui| about(ui, node, info, &f.labels[k], group));
        }
    }

    fn export_bar(&mut self, ui: &mut Ui, data: &Arc<NetworkData>, plan: &Arc<Plan>, folder: &str) {
        let save = |name: &str, filter: (&str, &str)| rfd::FileDialog::new().set_file_name(name).add_filter(filter.0, &[filter.1]).save_file();
        ui.horizontal_wrapped(|ui| {
            ui.label(crate::theme::bold("Export"));
            egui::ComboBox::from_id_salt("net_format").selected_text(self.format.label()).show_ui(ui, |ui| {
                for f in figure::Format::ALL {
                    ui.selectable_value(&mut self.format, f, f.label());
                }
            });
            if !self.format.is_vector() {
                ui.add(egui::DragValue::new(&mut self.dpi).range(72..=1200).suffix(" dpi"));
            }
            if ui.button("Save figure…").clicked() {
                let format = self.format;
                if let Some(path) = save(&format!("crosslink_network.{}", format.extension()), (format.label(), format.extension())) {
                    let svg = network::svg(data, &self.ann, plan, &self.style);
                    let opaque = self.style.background != Background::Transparent;
                    self.message = Some(match figure::export(&svg, format, self.dpi, opaque, &path) {
                        Ok(()) => (false, format!("Saved {}", path.display())),
                        Err(e) => (true, format!("Could not save {}: {e}", path.display())),
                    });
                }
            }
            if ui
                .button("Cytoscape…")
                .on_hover_text("The network for Cytoscape (.cyjs: File > Import > Network from File), with its layout and every attribute looked up")
                .clicked()
            {
                if let Some(path) = save("crosslink_network.cyjs", ("Cytoscape JSON", "cyjs")) {
                    self.message = Some(write(&path, network::cytoscape(data, &self.ann, plan, &self.style)));
                }
            }
            if ui.button("GraphML…").on_hover_text("For Cytoscape, Gephi, yEd, igraph or NetworkX").clicked() {
                if let Some(path) = save("crosslink_network.graphml", ("GraphML", "graphml")) {
                    self.message = Some(write(&path, network::graphml(data, &self.ann, plan, &self.style)));
                }
            }
            if ui.button("Tables…").on_hover_text("The proteins and the crosslinked pairs, as two tab-separated tables").clicked() {
                if let Some(dir) = rfd::FileDialog::new().pick_folder() {
                    let proteins = write(&dir.join("network_proteins.tsv"), network::node_table(data, &self.ann, plan, &self.style));
                    let pairs = write(&dir.join("network_pairs.tsv"), network::edge_table(data, &self.ann, plan, &self.style));
                    self.message = Some(if proteins.0 { proteins } else { pairs });
                }
            }
            let dir = Path::new(folder).join(network::FOLDER);
            if ui
                .add_enabled(!folder.is_empty(), egui::Button::new("All to results folder"))
                .on_hover_text(format!("SVG, PDF, PNG, this style, Cytoscape, GraphML and the tables into {}", dir.display()))
                .clicked()
            {
                self.message = Some(match network::write_all(data, &self.ann, plan, &self.style, &dir) {
                    Ok(files) => (false, format!("Wrote {} files to {}", files.len(), dir.display())),
                    Err(e) => (true, format!("Could not write the network: {e}")),
                });
            }
        });
        if let Some((error, text)) = &self.message {
            let color = if *error { Color32::from_rgb(210, 60, 60) } else { ui.visuals().weak_text_color() };
            ui.label(RichText::new(text).color(color));
        }
    }
}

fn write(path: &Path, text: String) -> (bool, String) {
    match std::fs::write(path, text) {
        Ok(()) => (false, format!("Saved {}", path.display())),
        Err(e) => (true, format!("Could not save {}: {e}", path.display())),
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use egui_kittest::Harness;

    /// Dragging a protein in the preview moves it there (and only it), on a
    /// real crosslink results folder (`SAGE_GUI_RESULTS`); the window after:
    /// `network_window_dragged.png` in `SAGE_GUI_FIGURES`
    #[test]
    #[ignore]
    fn network_window_drag() {
        let folder = std::env::var("SAGE_GUI_RESULTS").expect("set SAGE_GUI_RESULTS");
        let out = std::path::PathBuf::from(std::env::var("SAGE_GUI_FIGURES").expect("set SAGE_GUI_FIGURES"));
        let results = Arc::new(Results::load(Path::new(&folder)).unwrap());
        let view = NetworkView { open: true, ..Default::default() };
        let mut installed = false;
        let mut harness = Harness::builder().with_size(egui::vec2(1300.0, 1000.0)).build_state(
            move |ctx, view: &mut NetworkView| {
                if !installed {
                    crate::theme::install(ctx, false);
                    installed = true;
                    return;
                }
                egui::Window::new("Protein network")
                    .default_size([1260.0, 960.0])
                    .default_pos([10.0, 10.0])
                    .show(ctx, |ui| view.ui(ui, &results, &folder));
            },
            view,
        );
        harness.run_steps(6);
        // the most crosslinked protein, where it is drawn
        let (rect, f, accession, before, others) = {
            let v = harness.state();
            let (data, plan) = (v.data.as_ref().unwrap().0.clone(), v.plan.as_ref().unwrap().0.clone());
            let f = network::frame(&data, &v.ann, &plan, &v.style);
            let k = 0;
            let accession = data.nodes[plan.nodes[k]].accession.clone();
            let others: Vec<(f64, f64, f64)> = f.nodes.iter().enumerate().filter(|(i, _)| *i != k).map(|(_, n)| *n).collect();
            (v.preview_rect.unwrap(), f.clone(), accession, f.nodes[k], others)
        };
        let to_screen = |(x, y): (f64, f64)| {
            egui::pos2(rect.min.x + (x / f.width) as f32 * rect.width(), rect.min.y + (y / f.height) as f32 * rect.height())
        };
        let from = to_screen((before.0, before.1));
        let to = from + egui::vec2(-140.0, -90.0);
        let button = |pos, pressed| egui::Event::PointerButton {
            pos,
            button: egui::PointerButton::Primary,
            pressed,
            modifiers: egui::Modifiers::NONE,
        };
        harness.input_mut().events.push(egui::Event::PointerMoved(from));
        harness.run_steps(1);
        harness.input_mut().events.push(button(from, true));
        harness.run_steps(1);
        for i in 1..=10 {
            harness.input_mut().events.push(egui::Event::PointerMoved(from + (to - from) * (i as f32 / 10.0)));
            harness.run_steps(1);
        }
        harness.input_mut().events.push(button(to, false));
        harness.run_steps(3);
        let v = harness.state();
        assert!(v.dragging.is_none());
        assert_eq!(v.style.positions.keys().collect::<Vec<_>>(), [&accession], "only the protein dragged is moved");
        let (data, plan) = (v.data.as_ref().unwrap().0.clone(), v.plan.as_ref().unwrap().0.clone());
        let after = network::frame(&data, &v.ann, &plan, &v.style);
        let now = to_screen((after.nodes[0].0, after.nodes[0].1));
        assert!((now - to).length() < 3.0, "the protein is where it was dropped: {now:?} for {to:?}");
        let still: Vec<(f64, f64, f64)> = after.nodes.iter().skip(1).copied().collect();
        assert_eq!(still, others, "the others stay");
        std::fs::create_dir_all(&out).unwrap();
        let path = out.join("network_window_dragged.png");
        harness.render().expect("render").save(&path).unwrap();
        println!("{}", path.display());
    }

    /// The window on a real crosslink results folder (`SAGE_GUI_RESULTS`),
    /// with what was looked up for it if kept there: `network_window_*.png`
    /// in `SAGE_GUI_FIGURES`
    #[test]
    #[ignore]
    fn network_window_preview() {
        let folder = std::env::var("SAGE_GUI_RESULTS").expect("set SAGE_GUI_RESULTS");
        let out = std::path::PathBuf::from(std::env::var("SAGE_GUI_FIGURES").expect("set SAGE_GUI_FIGURES"));
        std::fs::create_dir_all(&out).unwrap();
        let results = Arc::new(Results::load(Path::new(&folder)).unwrap());
        let variants: [(&str, bool, Box<dyn Fn(&mut NetworkStyle)>); 4] = [
            ("complexes", false, Box::new(|_| {})),
            ("modules_dark", true, Box::new(|s| {
                s.group_by = GroupBy::Module;
                s.color_by = EdgeColor::Group;
                s.numbers = true;
            })),
            ("pathways", false, Box::new(|s| {
                s.group_by = GroupBy::Pathway;
                s.preferred = vec!["R-SCE-674695".into()];
                s.color_by = EdgeColor::Complex;
                s.layout = Layout::Groups;
                s.label_place = LabelPlace::Center;
                s.node_size = 22.0;
                s.node_size_max = 40.0;
            })),
            ("location_circle", false, Box::new(|s| {
                s.group_by = GroupBy::Location;
                s.layout = Layout::Circle;
                s.legend_side = LegendSide::Below;
                s.title = "RNA polymerase II".into();
            })),
        ];
        for (name, dark, change) in variants {
            let mut view = NetworkView { open: true, ..Default::default() };
            change(&mut view.style);
            let (results, folder) = (results.clone(), folder.clone());
            let mut installed = false;
            let mut harness = Harness::builder()
                .with_size(egui::vec2(1300.0, 1000.0))
                .build(move |ctx| {
                    if !installed {
                        crate::theme::install(ctx, dark);
                        crate::theme::set_dark(ctx, dark);
                        installed = true;
                        return;
                    }
                    egui::Window::new("Protein network")
                        .default_size([1260.0, 960.0])
                        .default_pos([10.0, 10.0])
                        .show(ctx, |ui| view.ui(ui, &results, &folder));
                });
            harness.run_steps(6);
            let path = out.join(format!("network_window_{name}.png"));
            harness.render().expect("render").save(&path).unwrap();
            println!("{}", path.display());
        }
    }
}
