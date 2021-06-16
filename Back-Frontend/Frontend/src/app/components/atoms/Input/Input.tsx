import React, { ChangeEvent, useState } from 'react';
import './Input.style.scss';

const Input = (props) => {
  const { label = '', type = 'text', value = '', onChange = () => {} } = props;
  const [inputValue, setInputValue] = useState(value);
  const handleChange = (event: ChangeEvent<HTMLInputElement>) => {
    const newVal = event.target?.value;
    setInputValue(newVal);
    onChange(newVal);
  };
  return (
    <div className="form-element input-container">
      <label className="form__label">{label}</label>
      <input
        className="form-input"
        type={type}
        value={inputValue}
        onChange={handleChange}
        placeholder={`Enter a ${label}...`}
      />
    </div>
  );
};

export default Input;
