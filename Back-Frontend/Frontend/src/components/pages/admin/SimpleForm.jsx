import React from 'react';
import { reduxForm, Field } from 'redux-form';
import { newField } from '../../common/newField';

const required = (value) => (value ? undefined : 'Required');

// const newField = ({
//   input,
//   placeholder,
//   type,
//   id,
//   meta: { touched, error },
//   ...rest
// }) => {
//   return (
//     <div>
//       <input {...input} type={type} placeholder={placeholder} id={id} />
//       {touched && error && <p style={{ color: "red" }}>{error}</p>}
//     </div>
//   );
// };

const SimpleForm = ({ handleSubmit, reset, pristine, submitting, valid }) => {
  return (
    <div className="col-sm-6">
      <div className="box box-info">
        <div className="box-header with-border">
          <h3 className="box-title">Redux Form</h3>
        </div>
        <form
          className="form-horizontal"
          onSubmit={handleSubmit((val) => console.log(val))}
        >
          <div className="box-body">
            <div className="form-group">
              <label htmlFor="first-name" className="col-sm-3 control-label">
                first name
              </label>
              <div className="col-sm-9">
                <Field
                  placeholder="enter first name"
                  type="text"
                  id="first-name"
                  name="first-name"
                  component={newField}
                  validate={[required]}
                />
              </div>
            </div>

            <button type="submit" disabled={!valid || pristine || submitting}>
              Submit
            </button>
            <button type="button" onClick={reset}>
              reset
            </button>
          </div>
        </form>
      </div>
    </div>
  );
};

export default reduxForm({ form: 'simpleForm' })(SimpleForm);
