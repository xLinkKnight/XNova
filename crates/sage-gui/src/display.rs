//! How the window is drawn, what happens when it cannot be, and a record of
//! unexpected stops.
//!
//! sage-gui draws with Direct3D 12 (through wgpu) by default. It needs no
//! OpenGL driver, which Remote Desktop sessions and virtual machines often
//! lack (`egui_glow requires opengl 2.0+`), and with no graphics card that
//! can draw the window it uses Windows' own software renderer (WARP), which
//! every Remote Desktop session has. OpenGL, which sage-gui used up to
//! v1.9.10, stays available: `XNOVA_RENDERER=opengl`.
//!
//! When the window cannot be opened, sage-gui starts again with the next way
//! of drawing (Direct3D 12, then the software renderer, then OpenGL), and
//! when none works it says why in a message box. A panic is written to
//! `%APPDATA%\sage-gui\crash.log`, and one that ends the window is also
//! shown, so that a window that closed can be explained.

use eframe::egui_wgpu;
use eframe::wgpu;
use std::sync::Arc;

/// The ways the window can be drawn, in the order they are tried
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Renderer {
    /// Direct3D 12 on a graphics card (else the software renderer)
    Direct3d,
    /// Direct3D 12 on Windows' software renderer (WARP)
    Software,
    /// OpenGL 2 or later (egui_glow)
    OpenGl,
}

/// Environment variable naming the renderer: `direct3d` (default),
/// `software` or `opengl`
pub const RENDERER_VAR: &str = "XNOVA_RENDERER";
/// Set when sage-gui started itself again with the next renderer
const FALLBACK_VAR: &str = "XNOVA_RENDERER_FALLBACK";

impl Renderer {
    pub fn key(self) -> &'static str {
        match self {
            Renderer::Direct3d => "direct3d",
            Renderer::Software => "software",
            Renderer::OpenGl => "opengl",
        }
    }

    pub fn describe(self) -> &'static str {
        match self {
            Renderer::Direct3d => "Direct3D 12",
            Renderer::Software => "Windows' software renderer (Direct3D 12 WARP)",
            Renderer::OpenGl => "OpenGL",
        }
    }

    /// The renderer `text` names (`dx12` and `d3d12` are Direct3D too)
    pub fn parse(text: &str) -> Option<Self> {
        match text.trim().to_ascii_lowercase().as_str() {
            "direct3d" | "d3d12" | "dx12" | "wgpu" => Some(Renderer::Direct3d),
            "software" | "warp" => Some(Renderer::Software),
            "opengl" | "gl" | "glow" => Some(Renderer::OpenGl),
            _ => None,
        }
    }

    /// The renderer to start with: `XNOVA_RENDERER`, else Direct3D 12
    pub fn chosen() -> Self {
        std::env::var(RENDERER_VAR).ok().and_then(|v| Self::parse(&v)).unwrap_or(Renderer::Direct3d)
    }

    /// The renderer to try when this one cannot open the window
    fn next(self) -> Option<Self> {
        match self {
            Renderer::Direct3d => Some(Renderer::Software),
            Renderer::Software => Some(Renderer::OpenGl),
            Renderer::OpenGl => None,
        }
    }

    /// eframe's options for this renderer
    pub fn apply(self, options: &mut eframe::NativeOptions) {
        match self {
            Renderer::OpenGl => options.renderer = eframe::Renderer::Glow,
            Renderer::Direct3d | Renderer::Software => {
                options.renderer = eframe::Renderer::Wgpu;
                options.wgpu_options = wgpu_configuration(self == Renderer::Software);
            }
        }
    }
}

/// Direct3D 12 only (Vulkan drivers have crashed opening a device on some
/// Windows computers), on the adapter [`pick_adapter`] chooses
fn wgpu_configuration(software: bool) -> egui_wgpu::WgpuConfiguration {
    egui_wgpu::WgpuConfiguration {
        wgpu_setup: egui_wgpu::WgpuSetup::CreateNew(egui_wgpu::WgpuSetupCreateNew {
            instance_descriptor: wgpu::InstanceDescriptor {
                backends: wgpu::Backends::DX12,
                flags: wgpu::InstanceFlags::from_build_config().with_env(),
                backend_options: wgpu::BackendOptions::from_env_or_default(),
            },
            power_preference: wgpu::PowerPreference::LowPower,
            native_adapter_selector: Some(Arc::new(move |adapters, surface| {
                let displays = display_adapters();
                pick_adapter(adapters, surface, software, &displays)
            })),
            device_descriptor: Arc::new(|adapter| wgpu::DeviceDescriptor {
                label: Some("XNova window"),
                required_features: wgpu::Features::empty(),
                // the largest pictures the adapter can hold (Circos and
                // network previews are as large as the window)
                required_limits: wgpu::Limits::default().using_resolution(adapter.limits()),
                memory_hints: wgpu::MemoryHints::MemoryUsage,
            }),
            trace_path: None,
        }),
        ..Default::default()
    }
}

/// The adapter to draw with: of those that can show the window, a graphics
/// card driving a display first (the window appears there), then any other
/// graphics card, then the software renderer; only the software renderer
/// when `software`
fn pick_adapter(
    adapters: &[wgpu::Adapter],
    surface: Option<&wgpu::Surface<'_>>,
    software: bool,
    displays: &[String],
) -> Result<wgpu::Adapter, String> {
    let rank = |a: &wgpu::Adapter| {
        let info = a.get_info();
        let cpu = info.device_type == wgpu::DeviceType::Cpu;
        let display = displays.iter().any(|d| d.eq_ignore_ascii_case(info.name.trim()));
        (cpu, !display)
    };
    let mut usable: Vec<&wgpu::Adapter> = adapters
        .iter()
        .filter(|a| surface.is_none_or(|s| a.is_surface_supported(s)))
        .filter(|a| !software || a.get_info().device_type == wgpu::DeviceType::Cpu)
        .collect();
    usable.sort_by_key(|a| rank(a));
    match usable.first() {
        Some(a) => Ok((*a).clone()),
        None => {
            let names: Vec<String> = adapters.iter().map(|a| a.get_info().name).collect();
            Err(format!(
                "no graphics adapter can draw the window{} (Direct3D 12 adapters: {})",
                if software { " with the software renderer" } else { "" },
                if names.is_empty() { "none".to_string() } else { names.join(", ") }
            ))
        }
    }
}

/// The names of the graphics cards with a display attached
#[cfg(windows)]
fn display_adapters() -> Vec<String> {
    use windows::Win32::Graphics::Dxgi::{CreateDXGIFactory1, IDXGIFactory1};
    let Ok(factory) = (unsafe { CreateDXGIFactory1::<IDXGIFactory1>() }) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for i in 0.. {
        let Ok(adapter) = (unsafe { factory.EnumAdapters1(i) }) else { break };
        if unsafe { adapter.EnumOutputs(0) }.is_err() {
            continue;
        }
        if let Ok(desc) = unsafe { adapter.GetDesc1() } {
            let len = desc.Description.iter().position(|&c| c == 0).unwrap_or(desc.Description.len());
            out.push(String::from_utf16_lossy(&desc.Description[..len]).trim().to_string());
        }
    }
    out
}

#[cfg(not(windows))]
fn display_adapters() -> Vec<String> {
    Vec::new()
}

/// The window could not be opened with `renderer`: start sage-gui again with
/// the next renderer (when the first came from sage-gui itself), else say
/// why. Returns true when another sage-gui was started.
pub fn window_failed(renderer: Renderer, error: &str) -> bool {
    let chosen_by_user = std::env::var_os(RENDERER_VAR).is_some() && std::env::var_os(FALLBACK_VAR).is_none();
    let mut tried: Vec<String> = std::env::var(FALLBACK_VAR)
        .ok()
        .map(|v| v.split(';').filter(|s| !s.is_empty()).map(str::to_string).collect())
        .unwrap_or_default();
    tried.push(format!("{}: {error}", renderer.describe()));
    if !chosen_by_user {
        if let Some(next) = renderer.next() {
            if let Ok(exe) = std::env::current_exe() {
                let started = std::process::Command::new(exe)
                    .args(std::env::args_os().skip(1))
                    .env(RENDERER_VAR, next.key())
                    .env(FALLBACK_VAR, tried.join(";"))
                    .spawn();
                if started.is_ok() {
                    return true;
                }
            }
        }
    }
    let advice = if chosen_by_user {
        format!(
            "{RENDERER_VAR} is set to {}. Remove it, or set it to direct3d, software or opengl, and start XNova again.",
            renderer.key()
        )
    } else {
        "Update the graphics driver. In a Remote Desktop session, the host's group policy \"Use hardware \
         graphics adapters for all Remote Desktop Services sessions\" lets programs use its graphics card."
            .to_string()
    };
    message_box(
        "XNova could not open its window",
        &format!("XNova could not open its window.\n\n{}\n\n{advice}", tried.join("\n")),
    );
    false
}

/// Record panics in `%APPDATA%\sage-gui\crash.log`; one that ends the window
/// (on the main thread) is also shown, so that a window that closed can be
/// explained. Background threads that catch their panics (reading a damaged
/// spectrum file) are only recorded.
pub fn record_panics() {
    let default = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |info| {
        let thread = std::thread::current();
        let name = thread.name().unwrap_or("unnamed").to_string();
        let text = format!(
            "{} XNova {} ({}), thread {name}: {info}\n{}\n",
            {
                let now = crate::clock::Stamp::now(std::time::Duration::ZERO);
                format!("{} {}", now.short_date(), now.time())
            },
            env!("CARGO_PKG_VERSION"),
            Renderer::chosen().describe(),
            std::backtrace::Backtrace::force_capture()
        );
        let path = crash_log();
        if let Some(path) = &path {
            if let Some(dir) = path.parent() {
                let _ = std::fs::create_dir_all(dir);
            }
            use std::io::Write;
            if let Ok(mut f) = std::fs::OpenOptions::new().create(true).append(true).open(path) {
                let _ = f.write_all(text.as_bytes());
            }
        }
        if name == "main" {
            let place = path.map(|p| format!("\n\nDetails are in {}.", p.display())).unwrap_or_default();
            message_box("XNova stopped", &format!("XNova stopped unexpectedly:\n\n{info}{place}"));
        }
        default(info);
    }));
}

/// `%APPDATA%\sage-gui\crash.log`
pub fn crash_log() -> Option<std::path::PathBuf> {
    let base = std::env::var_os("APPDATA").or_else(|| std::env::var_os("HOME"))?;
    Some(std::path::PathBuf::from(base).join("sage-gui").join("crash.log"))
}

#[cfg(windows)]
fn message_box(title: &str, text: &str) {
    use windows::core::PCWSTR;
    use windows::Win32::UI::WindowsAndMessaging::{MessageBoxW, MB_ICONERROR, MB_OK};
    let wide = |s: &str| s.encode_utf16().chain(std::iter::once(0)).collect::<Vec<u16>>();
    let (title, text) = (wide(title), wide(text));
    unsafe { MessageBoxW(None, PCWSTR(text.as_ptr()), PCWSTR(title.as_ptr()), MB_OK | MB_ICONERROR) };
}

#[cfg(not(windows))]
fn message_box(title: &str, text: &str) {
    eprintln!("{title}: {text}");
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn renderers_are_named_and_tried_in_turn() {
        assert_eq!(Renderer::parse(" OpenGL "), Some(Renderer::OpenGl));
        assert_eq!(Renderer::parse("dx12"), Some(Renderer::Direct3d));
        assert_eq!(Renderer::parse("warp"), Some(Renderer::Software));
        assert_eq!(Renderer::parse("vulkan"), None);
        assert_eq!(Renderer::Direct3d.next(), Some(Renderer::Software));
        assert_eq!(Renderer::Software.next(), Some(Renderer::OpenGl));
        assert_eq!(Renderer::OpenGl.next(), None);
        for r in [Renderer::Direct3d, Renderer::Software, Renderer::OpenGl] {
            assert_eq!(Renderer::parse(r.key()), Some(r));
        }
    }

    /// This computer's Direct3D 12 adapters: a graphics card is chosen over
    /// the software renderer, and the software renderer when asked for
    #[cfg(windows)]
    #[test]
    fn adapters_are_chosen_for_the_window() {
        let instance = wgpu::Instance::new(&wgpu::InstanceDescriptor { backends: wgpu::Backends::DX12, ..Default::default() });
        let adapters = instance.enumerate_adapters(wgpu::Backends::DX12);
        let names: Vec<String> = adapters.iter().map(|a| a.get_info().name).collect();
        eprintln!("Direct3D 12 adapters: {names:?}; with a display: {:?}", display_adapters());
        let Ok(software) = pick_adapter(&adapters, None, true, &[]) else {
            eprintln!("no software renderer here");
            return;
        };
        assert_eq!(software.get_info().device_type, wgpu::DeviceType::Cpu);
        let any = pick_adapter(&adapters, None, false, &display_adapters()).expect("an adapter");
        if adapters.iter().any(|a| a.get_info().device_type != wgpu::DeviceType::Cpu) {
            assert_ne!(any.get_info().device_type, wgpu::DeviceType::Cpu, "a graphics card first");
        }
        assert!(pick_adapter(&[], None, false, &[]).unwrap_err().contains("none"));
    }
}
