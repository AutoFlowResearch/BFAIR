import React from 'react';
import Select from 'react-select';
import AddOption from '../AddOption';
import './SelectBar.scss';

const SelectBar = (props) => {
  const {
    label,
    placeHolderText,
    selectedOption,
    options,
    onChange = () => {},
    multiple,
  } = props;
  const handleChange = (selectedOption) => {
    console.log(selectedOption);
  };
  return (
    <div className="form-element">
      <div className="label-container">
        <label className="form__label">{label}</label>
        <AddOption label={'Add New'} />
      </div>
      <Select
        isMulti={multiple}
        placeholder={`Enter${multiple ? '' : ' a'} ${placeHolderText}...`}
        onChange={handleChange}
        options={options}
      />
    </div>
  );
};

export default SelectBar;
