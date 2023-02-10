"use strict";

let molstar;

function init_results(structure_url, id) {
    (async () => {
        molstar = await MolstarPartialCharges.create("root");
        await load(structure_url, id);
    })().then(
        () => {
            // TODO: remove
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

    if (id === first_example) {
        document.getElementById("view_surface").setAttribute("checked", "true");
        await molstar.type.surface();
    } else {
        document
            .getElementById("colors_relative")
            .setAttribute("checked", "true");
        await molstar.type.default();
    }

    updateRelativeColor();
    mountTypeControls();
    mountColorControls();
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
    const relative = document.getElementById("colors_relative");
    const absolute = document.getElementById("colors_absolute");
    const range = document.getElementById("max_value");
    if (!structure || !relative || !absolute) return;
    structure.onclick = async () => await updateDefaultColor();
    relative.onclick = async () => await updateRelativeColor();
    absolute.onclick = async () => await updateAbsoluteColor();
    range.oninput = async () => await updateRange();
}

async function updateDefaultColor() {
    const input = document.getElementById("max_value");
    if (!input) return;
    input.setAttribute("disabled", "true");
    await molstar.color.default();
}

async function updateRelativeColor() {
    const input = document.getElementById("max_value");
    if (!input) return;
    input.setAttribute("disabled", "true");
    const charge = await molstar.charges.getRelativeCharge();
    input.value = charge.toFixed(3);
    // ? updates viewer each time the radio button is clicked ?
    await molstar.color.relative();
}

async function updateAbsoluteColor() {
    const input = document.getElementById("max_value");
    if (!input) return;
    input.removeAttribute("disabled");
    await molstar.color.relative();
}

async function updateRange() {
    const input = document.getElementById("max_value");
    if (!input) return;
    const value = Number(input.value);
    if (isNaN(value)) return;
    await molstar.color.absolute(value);
}
