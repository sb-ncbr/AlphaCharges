'use strict';


const spinner = '<span class="spinner-border spinner-border-sm" role="status" ' +
    'aria-hidden="true" style="animation-duration: 1.5s"></span>';

function update_litemol_colors(min_color, max_color) {
    LiteMolChargesViewerEventQueue.send("lm-set-default-color-scheme", {
        minVal: min_color,
        maxVal: max_color,
        fallbackColor: {r: 0, g: 255, b: 0},
        minColor: {r: 255, g: 0, b: 0},
        maxColor: {r: 0, g: 0, b: 255},
        middleColor: {r: 255, g: 255, b: 255}
    });
}


function init_results(idd, chg_range) {

    const $select = idd

    let $min_value = $('#min_value');
    let $max_value = $('#max_value');

        const id = $select;
        $.ajax({
            url: get_format_url,
            success: function (format) {
                LiteMolChargesViewerEventQueue.send("lm-load-molecule", {
                    structure_url: get_structure_url,
                    charges_url: get_charges_url,
                    structure_format: format,
                    charges_format: 'TXT'
                });
            }
        })




        $('input:radio[name=colors]').prop('disabled', false);

        if ($('input[name=colors]:checked').val() === 'Relative') {
            $min_value.val(-chg_range);
            $max_value.val(chg_range);
            $min_value.trigger('input');
        }


    $('#min_value, #max_value').on('input', function () {
        update_litemol_colors(parseFloat($('#min_value').val()), parseFloat($('#max_value').val()));
        $min_value.attr('max', $max_value.val());
        $max_value.attr('min', $min_value.val());
    });

    let $colors = $('input[name=colors]');
    $colors.on('change', function () {
        let coloring = $('input[name=colors]:checked').val();
        if (coloring === 'Relative') {
            LiteMolChargesViewerEventQueue.send('lm-use-default-themes', {value: false});
            const id = $select;
            $min_value.val(-chg_range);
            $max_value.val(chg_range);

            update_litemol_colors(null, null);
            $min_value.prop('disabled', true);
            $max_value.prop('disabled', true);
        } else if (coloring === 'Absolute') {
            LiteMolChargesViewerEventQueue.send('lm-use-default-themes', {value: false});
            $min_value.prop('disabled', false);
            $max_value.prop('disabled', false);
            $min_value.trigger('input');
        } else {
            /* Coloring by elements */
            LiteMolChargesViewerEventQueue.send('lm-use-default-themes', {value: true});

        }
    });


    let $view = $('input[name=view]');
    $view.on('change', function () {
        let v = $('input[name=view]:checked').val();
        if (v === 'Cartoon') {
            LiteMolChargesViewerEventQueue.send('lm-switch-to-cartoons');
        } else if (v === 'Balls and sticks') {
            LiteMolChargesViewerEventQueue.send('lm-switch-to-bas');
        } else {
            /* Surface */
            LiteMolChargesViewerEventQueue.send('lm-switch-to-surface');
        }
    });
    $colors.filter(':checked').trigger('change');
}


