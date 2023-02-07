"use strict";

let charge;
let molstar;

function init_results(structure_url, id) {
    (async () => {
        molstar = await MolstarPartialCharges.create("root");
        await load(structure_url, id);
    })().then(
        () => {
            console.log("Mol* initialization ✅");
        },
        (error) => {
            console.log("Mol* initialization ❌", error);
        }
    );
}

async function load(structure_url, id) {
    const first_example = "P34712_7.2_4";

    await molstar.load(structure_url);

    const cartoonOff = switchOffCartoonView();
    if (cartoonOff) {
        await molstar.type.ballAndStick();
        document.getElementById("view_bas").setAttribute("checked", "true");
    } else if (id === first_example) {
        await molstar.type.surface();
        document.getElementById("view_surface").setAttribute("checked", "true");
    } else {
        await molstar.type.default();
        document.getElementById("view_cartoon").setAttribute("checked", "true");
    }
    document.getElementById("colors_charges").setAttribute("checked", "true");

    let maxAbsoluteRelativeCharge = Number(
        molstar.charges.getRelativeCharge().toFixed(3)
    );
    await updateSliderMax(maxAbsoluteRelativeCharge);
    await updateCharge(maxAbsoluteRelativeCharge);

    mountTypeControls();
    mountColorControls();
    mountRangeControls();
    addEventListeners();
}

function switchOffCartoonView() {
    const view = document.getElementById("view_cartoon");
    if (!view) return false;
    if (!molstar.type.isDefaultApplicable()) {
        view.setAttribute("disabled", "true");
        return true;
    } else {
        view.removeAttribute("disabled");
        return false;
    }
}

async function updateSliderMax(max) {
    const slider = document.getElementById("input_range_max_charge");
    if (!slider) return;
    const old_value = Number(slider.value);
    slider.max = max.toFixed(3);
    if (old_value > max) {
        await updateCharge(max);
    }
}

async function updateCharge(chg) {
    if (isNaN(chg)) return;
    const header = document.getElementById("input_label_max_charge");
    const slider = document.getElementById("input_range_max_charge");
    if (!header || !slider) return;
    charge = Number(chg.toFixed(3));
    header.innerText = `Charge range: ${charge.toFixed(3)}`;
    slider.value = `${charge.toFixed(3)}`;
    await molstar.color.absolute(charge);
}

function mountTypeControls() {
    const cartoon = document.getElementById("view_cartoon");
    const surface = document.getElementById("view_surface");
    const bas = document.getElementById("view_bas");
    if (!cartoon || !surface || !bas) return;
    cartoon.onclick = async () => await molstar.type.default();
    surface.onclick = async () => await molstar.type.surface();
    bas.onclick = async () => await molstar.type.ballAndStick();
}

function mountColorControls() {
    const structure = document.getElementById("colors_structure");
    const absolute = document.getElementById("colors_charges");
    if (!structure || !absolute) return;
    structure.onclick = async () => await molstar.color.default();
    absolute.onclick = async () => {
        const relativeCharge = molstar.charges.getRelativeCharge();
        await updateSliderMax(relativeCharge);
        await updateCharge(relativeCharge);
    };
}

function mountRangeControls() {
    const slider = document.getElementById("input_range_max_charge");
    const input = document.getElementById("input_text_max_charge");
    if (!slider || !input) return;
    slider.oninput = async () => {
        const slider = document.getElementById("input_range_max_charge");
        if (!slider) return;
        const chg = Number(slider.value);
        await updateCharge(chg);
    };
    input.oninput = async () => {
        const input = document.getElementById("input_text_max_charge");
        if (!input) return;
        const value = input.value;
        const placeholder = input.placeholder;
        const defaultValue = !isNaN(Number(placeholder))
            ? Number(placeholder)
            : 0;
        const parsed =
            value !== "" && !isNaN(Number(value))
                ? Number(value)
                : defaultValue;
        await updateSliderMax(parsed);
    };
}
